module VtkReader

using Base64: base64decode
using CodecZlib: ZlibDecompressor
using LightXML:
    LightXML, XMLElement, parse_string, attribute, child_elements, free, find_element

export read_vtk

const HEADER_TYPE = UInt64
const DATATYPE_MAPPING = Dict("Float64" => Float64)

struct DataArray{T}
    type::DataType
    name::String
    number_of_components::Int
    number_of_tuples::Int
    format::Symbol
    offset::Int
    function DataArray(xml::XMLElement)
        LightXML.name(xml) !== "DataArray" && error("Need a DataArray!")
        type = DATATYPE_MAPPING[attribute(xml, "type"; required=true)]
        name = attribute(xml, "Name"; required=true)
        noc = parse(Int, attribute(xml, "NumberOfComponents"; required=true))
        _not = attribute(xml, "NumberOfTuples")
        not = isnothing(_not) ? 0 : parse(Int, _not)
        format = Symbol(attribute(xml, "format"; required=true))
        if format !== :appended
            error("Only appended data format valid!\n")
        end
        offset = parse(Int, attribute(xml, "offset"; required=true))
        new{type}(type, name, noc, not, format, offset)
    end
end

struct DataArrayMapping
    points::DataArray
    point_data::Vector{DataArray}
    field_data::Vector{DataArray}
    data::Vector{UInt8}
end

function get_xml_and_data(file::String)
    raw_file = read(file, String)

    # find appended data
    marker = findfirst("<AppendedData encoding=\"raw\">", raw_file)
    if isnothing(marker)
        error(
            "Invalid VTK-file!\n",
            "Can only read files with raw encoded data, specified by block:\n",
            "  <AppendedData encoding=\"raw\">\n    ...raw data...\n  </AppendedData>\n",
            "Could not find `<AppendedData encoding=\"raw\">` in file!\n",
        )
    end
    offset_begin_marker = findnext("_", raw_file, last(marker))
    if isnothing(offset_begin_marker)
        error(
            "Could not find the begin of the appended data!\n",
            "Usually this is marked after the _ character, which could not ",
            "be found after the <AppendedData encoding=\"raw\"> statement.\n",
        )
    end
    offset_begin = first(offset_begin_marker) + 1
    offset_end_marker = findnext("\n  </AppendedData>", raw_file, offset_begin)
    if isnothing(offset_end_marker)
        error(
            "Invalid VTK-file!\n",
            "Can only read files with raw encoded data, specified by block:\n",
            "  <AppendedData encoding=\"raw\">\n    ...raw data...\n  </AppendedData>\n",
            "Could not find `\\n  </AppendedData>` in file!\n",
        )
    end
    offset_end = first(offset_end_marker) - 1
    data = Vector{UInt8}(raw_file[offset_begin:offset_end])

    # xml contents
    xml_string = raw_file[begin:(offset_begin-1)] * raw_file[(offset_end+1):end]

    return xml_string, data
end

function extract_data_arrays(xml_string::String)

    # parse xml document
    xml_doc = parse_string(xml_string)
    root = LightXML.root(xml_doc)
    LightXML.name(root) !== "VTKFile" && error("File should be a valid VTKFile!")

    # get elements...
    unstructured_grid = find_element(root, "UnstructuredGrid")
    isnothing(unstructured_grid) && error("VTK-file has to contain a `UnstructuredGrid`!\n")

    piece = find_element(unstructured_grid, "Piece")
    isnothing(piece) && error("UnstructuredGrid does not contain `Piece`!\n")

    points = find_element(piece, "Points")
    isnothing(points) && error("Piece does not contain `Points`!\n")
    points_da_xml = find_element(points, "DataArray")
    isnothing(points_da_xml) && error("No DataArray found for `Points`!\n")
    points_da = DataArray(points_da_xml)

    point_data = find_element(piece, "PointData")
    isnothing(point_data) && error("Piece does not contain `PointData`!\n")
    point_data_arrays = Vector{DataArray}()
    for xml_elem in child_elements(point_data)
        push!(point_data_arrays, DataArray(xml_elem))
    end

    field_data = find_element(unstructured_grid, "FieldData")
    isnothing(field_data) && error("UnstructuredGrid does not contain `FieldData`!\n")
    field_data_arrays = Vector{DataArray}()
    for xml_elem in child_elements(field_data)
        push!(field_data_arrays, DataArray(xml_elem))
    end

    free(xml_doc)

    return points_da, point_data_arrays, field_data_arrays
end

function get_data_array_mapping(file::String)
    xml_string, data = get_xml_and_data(file)
    points, point_data, field_data = extract_data_arrays(xml_string)
    return DataArrayMapping(points, point_data, field_data, data)
end

function get_data(da::DataArray, data::Vector{UInt8})::VecOrMat{Float64}
    # extract number of bytes from header
    start = da.offset + 1
    stop = da.offset + 4 * sizeof(HEADER_TYPE)
    header = Int.(reinterpret(HEADER_TYPE, data[start:stop]))
    n_bytes = header[4]

    # get start and stop index for element in data
    start = da.offset + 4 * sizeof(HEADER_TYPE) + 1
    stop = start + n_bytes - 1

    # get the array out of compressed data
    data_decompressed = transcode(ZlibDecompressor, data[start:stop])
    data_float::Vector{da.type} = reinterpret(Float64, data_decompressed)

    if da.number_of_components > 1
        return reshape(data_float, da.number_of_components, :)
    else
        return data_float
    end
end

function get_dict(dam::DataArrayMapping)
    d = Dict{String, VecOrMat{Float64}}()
    d["Position"] = get_data(dam.points, dam.data)
    for pd in dam.point_data
        d[pd.name] = get_data(pd, dam.data)
    end
    for pd in dam.field_data
        d[pd.name] = get_data(pd, dam.data)
    end

    return d
end

"""
    read_vtk(file::String)

Read .vtu-file containing simulation results of a time step.

# Arguments
- `file::String`: path to VTK .vtu-file

# Returns
- `Dict{String, VecOrMat{Float64}}`: simulation results as a dictionary

# Examples

```julia-repl
julia> read_vtk("ExampleSimulation_t3000.vtu")
Dict{String, VecOrMat{Float64}} with 6 entries:
  "Position"     => [-0.497302 -0.497303 … 0.497303 0.497302; -0.0225023 -0…
  "Damage"       => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0…
  "Time"         => [0.000237879]
  "Displacement" => [0.00019766 0.000196793 … -0.000196793 -0.00019766; -2.…
  "Velocity"     => [-2.55436 -1.5897 … 1.5897 2.55436; 0.827107 0.234996 ……
  "ForceDensity" => [-7.00131e9 -8.45411e9 … 8.45411e9 7.00131e9; -5.0164e8…
```
"""
function read_vtk(file::String)
    if !endswith(file, ".vtu")
        _, extension = splitext(basename(file))
        msg = "cannot read file with extension $extension, specify a valid .vtu-file!"
        throw(AssertionError(msg))
    end
    dam = get_data_array_mapping(file)
    d = get_dict(dam)
    return d
end

end # module VtkReader
