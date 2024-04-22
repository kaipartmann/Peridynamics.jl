module VtkReader

using Base64: base64decode
using CodecZlib: ZlibDecompressor
using LightXML: LightXML, XMLDocument, XMLElement, parse_string, attribute, child_elements,
                free, find_element, get_elements_by_tagname

export read_vtk

const HEADER_TYPE = UInt64
const DATATYPE_MAPPING = Dict("Float64" => Float64, "Int64" => Int64)

abstract type AbstractDataArray end

struct DataArray{T,N} <: AbstractDataArray
    name::Symbol
    offset::Int
    function DataArray(xml::XMLElement)
        if LightXML.name(xml) != "DataArray"
            throw(ArgumentError("xml element has to be a DataArray!"))
        end
        type_str = attribute(xml, "type"; required=true)
        name_str = attribute(xml, "Name"; required=true)
        noc_str = attribute(xml, "NumberOfComponents"; required=true)
        offset_str = attribute(xml, "offset"; required=true)
        format_str = attribute(xml, "format"; required=true)
        if format_str != "appended"
            throw(ArgumentError("only DataArrays with format appended are valid!\n"))
        end
        if !haskey(DATATYPE_MAPPING, type_str)
            msg = "data of type $type_str for data $name_str not supported!\n"
            throw(ArgumentError(msg))
        end
        type = DATATYPE_MAPPING[type_str]
        name = Symbol(name_str)
        noc = parse(Int, noc_str)
        offset = parse(Int, offset_str)
        return new{type,noc}(name, offset)
    end
end

struct PDataArray{T,N} <: AbstractDataArray
    name::Symbol
    das::Vector{DataArray{T,N}}
    function PDataArray(das::Vector{DataArray{T,N}}) where {T,N}
        name = first(das).name
        return new{T,N}(name, das)
    end
end

element_type(::DataArray{T,N}) where {T,N} = T
element_type(::PDataArray{T,N}) where {T,N} = T

function get_raw_xml_and_data(file::String)
    raw_file = read(file, String)

    # find appended data
    marker = findfirst("<AppendedData encoding=\"raw\">", raw_file)
    if isnothing(marker)
        error("Invalid VTK-file!\n",
              "Can only read files with raw encoded data, specified by block:\n",
              "  <AppendedData encoding=\"raw\">\n    ...raw data...\n  </AppendedData>\n",
              "Could not find `<AppendedData encoding=\"raw\">` in file!\n")
    end
    offset_begin_marker = findnext("_", raw_file, last(marker))
    if isnothing(offset_begin_marker)
        error("Could not find the begin of the appended data!\n",
              "Usually this is marked after the _ character, which could not ",
              "be found after the <AppendedData encoding=\"raw\"> statement.\n")
    end
    offset_begin = first(offset_begin_marker) + 1
    offset_end_marker = findnext("\n  </AppendedData>", raw_file, offset_begin)
    if isnothing(offset_end_marker)
        error("Invalid VTK-file!\n",
              "Can only read files with raw encoded data, specified by block:\n",
              "  <AppendedData encoding=\"raw\">\n    ...raw data...\n  </AppendedData>\n",
              "Could not find `\\n  </AppendedData>` in file!\n")
    end
    offset_end = first(offset_end_marker) - 1
    raw_data = Vector{UInt8}(raw_file[offset_begin:offset_end])

    # xml contents
    raw_xml_str = raw_file[begin:(offset_begin-1)] * raw_file[(offset_end+1):end]

    return raw_xml_str, raw_data
end

function get_raw_xml_and_data(files::Vector{String})
    vec_raw_xml_str = Vector{String}(undef, length(files))
    vec_raw_data = Vector{Vector{UInt8}}(undef, length(files))
    Threads.@threads for i in eachindex(files)
        raw_xml_str, raw_data = get_raw_xml_and_data(files[i])
        vec_raw_xml_str[i] = raw_xml_str
        vec_raw_data[i] = raw_data
    end
    return vec_raw_xml_str, vec_raw_data
end

function find_data_xml(xml_doc::XMLDocument)
    root = LightXML.root(xml_doc)
    LightXML.name(root) !== "VTKFile" && error("File should be a valid VTKFile!")
    unstruct_grid_xml = find_element(root, "UnstructuredGrid")
    isnothing(unstruct_grid_xml) && error("VTK-file has to contain a `UnstructuredGrid`!\n")
    piece_xml = find_element(unstruct_grid_xml, "Piece")
    isnothing(piece_xml) && error("UnstructuredGrid does not contain `Piece`!\n")
    points_xml = find_element(piece_xml, "Points")
    isnothing(points_xml) && error("Piece does not contain `Points`!\n")
    position_xml = find_element(points_xml, "DataArray")
    isnothing(position_xml) && error("No DataArray found for `Points`!\n")
    point_data_xml = find_element(piece_xml, "PointData")
    isnothing(point_data_xml) && error("Piece does not contain `PointData`!\n")
    field_data_xml = find_element(unstruct_grid_xml, "FieldData")
    isnothing(field_data_xml) && error("UnstructuredGrid does not contain `FieldData`!\n")
    return position_xml, point_data_xml, field_data_xml
end

@inline extract_das(data_xml::XMLElement) = DataArray.(child_elements(data_xml))

function parse_data_arrays(raw_xml_str::AbstractString)
    xml_doc = parse_string(raw_xml_str)
    position_xml, point_data_xml, field_data_xml = find_data_xml(xml_doc)
    position_da = DataArray(position_xml)
    point_das = extract_das(point_data_xml)
    field_das = extract_das(field_data_xml)
    free(xml_doc)
    return position_da, point_das, field_das
end

function parse_data_arrays(vec_raw_xml_str::Vector{S}) where {S<:AbstractString}
    n = length(vec_raw_xml_str)
    vec_position_da = Vector{DataArray}(undef, n)
    vec_point_das = Vector{Vector{DataArray}}(undef, n)
    vec_n_point_das = zeros(Int, n)
    vec_field_das = Vector{Vector{DataArray}}(undef, n)
    vec_n_field_das = zeros(Int, n)
    Threads.@threads for i in eachindex(vec_raw_xml_str)
        position_da, point_das, field_das = parse_data_arrays(vec_raw_xml_str[i])
        vec_position_da[i] = position_da
        vec_point_das[i] = point_das
        vec_n_point_das[i] = length(point_das)
        vec_field_das[i] = field_das
        vec_n_field_das[i] = length(field_das)
    end
    position_pda = PDataArray(promote_type.(vec_position_da))

    allequal(vec_n_point_das) || error("corrupted pvtu file!\n")
    n_point_das = first(vec_n_point_das)
    point_pdas = Vector{PDataArray}(undef, n_point_das)
    Threads.@threads for i in 1:n_point_das
        point_pdas[i] = PDataArray([pdas[i] for pdas in vec_point_das])
    end

    allequal(vec_n_field_das) || error("corrupted pvtu file!\n")
    n_field_das = first(vec_n_field_das)
    field_pdas = Vector{PDataArray}(undef, n_field_das)
    Threads.@threads for i in 1:n_field_das
        field_pdas[i] = PDataArray([fdas[i] for fdas in vec_field_das])
    end

    return position_pda, point_pdas, field_pdas
end

function get_data(da::DataArray{T,N}, raw_data::Vector{UInt8}) where {T,N}
    data::Matrix{T} = reshape(reinterpret(T, get_decompressed_data(da, raw_data)), N, :)
    return data
end

function get_data(da::DataArray{T,1}, raw_data::Vector{UInt8}) where {T}
    data::Vector{T} = reinterpret(T, get_decompressed_data(da, raw_data))
    return data
end

function get_data(pda::PDataArray{T,N}, raw_data::Vector{Vector{UInt8}}) where {T,N}
    vec_data = get_data.(pda.das, raw_data)
    return reduce(hcat, vec_data)
end

function get_data(pda::PDataArray{T,1}, raw_data::Vector{Vector{UInt8}}) where {T}
    vec_data = get_data.(pda.das, raw_data)
    return reduce(vcat, vec_data)
end

function get_field_data(da::DataArray{T,N}, raw_data::Vector{UInt8}) where {T,N}
    return get_data(da, raw_data)
end

function get_field_data(pda::PDataArray{T,N}, raw_data::Vector{Vector{UInt8}}) where {T,N}
    return unique(get_data(pda, raw_data); dims=2)
end

function get_field_data(pda::PDataArray{T,1}, raw_data::Vector{Vector{UInt8}}) where {T}
    return unique(get_data(pda, raw_data))
end

function get_decompressed_data(da::DataArray{T,N}, raw_data::Vector{UInt8}) where {T,N}
    # extract number of bytes from header
    start = da.offset + 1
    stop = da.offset + 4 * sizeof(HEADER_TYPE)
    header = Int.(reinterpret(HEADER_TYPE, raw_data[start:stop]))
    n_bytes = header[4]

    # get start and stop index for element in data
    start = da.offset + 4 * sizeof(HEADER_TYPE) + 1
    stop = start + n_bytes - 1

    # get the array out of compressed data
    data_decompressed = transcode(ZlibDecompressor, raw_data[start:stop])

    return data_decompressed
end

function promote_vec_or_mat_elemtype(pda::D1, pdas::Vector{D2},
                                     fdas::Vector{D3}) where {D1,D2,D3<:AbstractDataArray}
    element_type_pda = element_type(pda)
    element_types_pdas = element_type.(pdas)
    element_types_fdas = element_type.(fdas)
    return typejoin(element_type_pda, element_types_pdas..., element_types_fdas...)
end

function get_dict(pda, pdas, fdas, raw_data)
    type = promote_vec_or_mat_elemtype(pda, pdas, fdas)
    if isconcretetype(type)
        d = Dict{Symbol,VecOrMat{type}}()
    else
        d = Dict{Symbol,VecOrMat{<:type}}()
    end
    d[:position] = get_data(pda, raw_data)
    for da in pdas
        d[da.name] = get_data(da, raw_data)
    end
    for da in fdas
        d[da.name] = get_field_data(da, raw_data)
    end
    return d
end

function get_vtu_files(file::AbstractString)
    xml_string = read(file, String)
    xml_doc = parse_string(xml_string)
    root = LightXML.root(xml_doc)
    LightXML.name(root) != "VTKFile" && error("File should be a valid VTKFile!")
    unstruct_grid = find_element(root, "PUnstructuredGrid")
    isnothing(unstruct_grid) && error("VTK-file has to contain a `PUnstructuredGrid`!\n")
    pieces = get_elements_by_tagname(unstruct_grid, "Piece")
    isnothing(pieces) && error("UnstructuredGrid does not contain `Piece`!\n")
    vtu_files::Vector{String} = attribute.(pieces, "Source"; required=true)
    add_pwd!(vtu_files, dirname(file))
    return vtu_files
end

function add_pwd!(vtu_files::Vector{String}, dirname_file::String)
    for i in eachindex(vtu_files)
        vtu_files[i] = joinpath(dirname_file, vtu_files[i])
    end
    return nothing
end

function _read_vtu(vtu_file::AbstractString)
    raw_xml_str, raw_data = get_raw_xml_and_data(vtu_file)
    position_da, point_das, field_das = parse_data_arrays(raw_xml_str)
    data_dict = get_dict(position_da, point_das, field_das, raw_data)
    return data_dict
end

function _read_pvtu(file::AbstractString)
    vtu_files = get_vtu_files(file)
    vec_raw_xml_str, vec_raw_data = get_raw_xml_and_data(vtu_files)
    position_pda, point_pdas, field_pdas = parse_data_arrays(vec_raw_xml_str)
    data_dict = get_dict(position_pda, point_pdas, field_pdas, vec_raw_data)
    return data_dict
end

is_vtu(file::AbstractString) = endswith(file, ".vtu") ? true : false
is_pvtu(file::AbstractString) = endswith(file, ".pvtu") ? true : false

"""
    read_vtk(file::AbstractString)

Read vtu or pvtu file containing simulation results of a time step.

# Arguments
- `file::String`: Path to VTK file in vtu or pvtu format

# Returns
- `Dict{String, VecOrMat{Float64}}`: Simulation results as a dictionary

# Examples

```julia-repl
julia> read_vtk("ExampleSimulation_t3000.vtu")
Dict{String, VecOrMat{Float64}} with 6 entries:
  "position"     => [-0.497302 -0.497303 … 0.497303 0.497302; -0.0225023 -0…
  "damage"       => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0…
  "time"         => [0.000237879]
  "displacement" => [0.00019766 0.000196793 … -0.000196793 -0.00019766; -2.…
  "velocity"     => [-2.55436 -1.5897 … 1.5897 2.55436; 0.827107 0.234996 ……
  "b_int"        => [-7.00131e9 -8.45411e9 … 8.45411e9 7.00131e9; -5.0164e8…
```
"""
function read_vtk(file::AbstractString)
    isfile(file) || throw(ArgumentError("cannot find file $(file)!\n"))
    if is_vtu(file)
        d = _read_vtu(file)
    elseif is_pvtu(file)
        d = _read_pvtu(file)
    else
        _, extension = splitext(basename(file))
        msg = "cannot read file with extension $extension, specify a valid file!\n"
        throw(ArgumentError(msg))
    end
    return d
end

end # module VtkReader
