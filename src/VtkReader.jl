module VtkReader

using Base64: base64decode
using CodecZlib: ZlibDecompressor
using LightXML:
    LightXML, XMLElement, parse_string, attribute, child_elements, free, find_element

export SimResult, read_vtk

const HEADER_TYPE = UInt64

struct SimResult
    position::Matrix{Float64}
    time::Float64
    damage::Vector{Float64}
    displacement::Matrix{Float64}
end

function Base.show(io::IO, ::MIME"text/plain", sr::SimResult)
    print(io, typeof(sr), " with fields:")
    for field in fieldnames(typeof(sr))
        field_val = getfield(sr, field)
        field_type = typeof(field_val)
        if !isempty(field_val)
            if field_type <: AbstractArray
                print(io, "\n  ", rpad(string(field) * ":", 14))
                Base.array_summary(io, field_val, axes(field_val))
            else
                print(io, "\n  ", rpad(string(field) * ":", 14), field_type)
            end
        end
    end
    return nothing
end

function get_xml_and_data(file::String)
    ## function read_vtk(file::String)
    raw_file = read(file, String)

    ## find appended data
    marker = findfirst("<AppendedData encoding=\"raw\">", raw_file)
    offset_begin = first(findnext("_", raw_file, last(marker))) + 1
    offset_end = first(findnext("</AppendedData>", raw_file, offset_begin)) - 1
    data = Vector{UInt8}(rstrip(raw_file[offset_begin:offset_end]))

    ## xml contents
    xml_contents = raw_file[1:(offset_begin - 1)] * "\n  </AppendedData>\n</VTKFile>"

    ## open end extract xml document
    xml_doc = parse_string(xml_contents)

    return xml_doc, data
end

function get_data_arrays(xml_doc)
    root = LightXML.root(xml_doc)
    @assert LightXML.name(root) == "VTKFile"
    points = root["UnstructuredGrid"][1]["Piece"][1]["Points"][1]
    point_data = root["UnstructuredGrid"][1]["Piece"][1]["PointData"][1]
    field_data = root["UnstructuredGrid"][1]["FieldData"][1]

    ## extract points
    position_da = find_element(points, "DataArray")

    ## extract time
    time_da = find_element(field_data, "DataArray")
    @assert LightXML.attribute(time_da, "Name"; required=true) == "time"

    ## get point data
    point_da_names = Vector{String}()
    point_da = Vector{XMLElement}()
    for xml_element in child_elements(point_data)
        @assert LightXML.name(xml_element) == "DataArray"
        push!(point_da_names, attribute(xml_element, "Name"; required=true))
        push!(point_da, xml_element)
    end

    return position_da, time_da, point_da_names, point_da
end

function get_data(xml_element, data)
    # Ensure the correct type of of the XML element
    @assert LightXML.name(xml_element) == "DataArray"

    # extract number of bytes from header
    offset = parse(Int, attribute(xml_element, "offset"; required=true))
    start = offset + 1
    stop = offset + 4 * sizeof(HEADER_TYPE)
    header = Int.(reinterpret(HEADER_TYPE, data[start:stop]))
    n_bytes = header[4]

    # get start and stop index for element in data
    start = offset + 4 * sizeof(HEADER_TYPE) + 1
    stop = start + n_bytes - 1

    # get the array out of compressed data
    data_decompressed = transcode(ZlibDecompressor, data[start:stop])
    data_float::Vector{Float64} = reinterpret(Float64, data_decompressed)

    return data_float
end

get_position(da, data) = reshape(get_data(da, data), 3, :)

get_time(da, data) = first(get_data(da, data))

function get_damage(pdan, pda, data)
    id = findfirst(x -> x == "damage", pdan)
    return get_data(pda[id], data)
end

function get_displacement(pdan, pda, data)
    id = findfirst(x -> x == "displacement", pdan)
    if !isnothing(id)
        displacement = reshape(get_data(pda[id], data), 3, :)
    else
        displacement = Array{Float64,2}(undef, 0, 0)
    end
    return displacement
end

function get_result(position_da, time_da, point_da_names, point_da, data)
    position = get_position(position_da, data)
    time = get_time(time_da, data)
    damage = get_damage(point_da_names, point_da, data)
    displacement = get_displacement(point_da_names, point_da, data)
    return SimResult(position, time, damage, displacement)
end

function read_vtk(file::String)
    if !endswith(file, ".vtu")
        _, extension = splitext(basename(file))
        msg = "cannot read file with extension $extension, specify a valid .vtu-file!"
        throw(AssertionError(msg))
    end
    xml_file, data = get_xml_and_data(file)
    position_da, time_da, point_da_names, point_da = get_data_arrays(xml_file)
    result = get_result(position_da, time_da, point_da_names, point_da, data)
    free(xml_file)
    return result
end

end # module VtkReader
