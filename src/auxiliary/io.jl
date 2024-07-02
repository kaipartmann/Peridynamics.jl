@inline default_export_fields() = [:displacement, :damage]

@inline function default_fields_spec(ms::AbstractMultibodySetup)
    fields_spec = Dict{Symbol,Vector{Symbol}}()
    for body_name in each_body_name(ms)
        fields_spec[body_name] = default_export_fields()
    end
    return fields_spec
end

function get_export_fields(body::AbstractBody, solver::AbstractTimeSolver,
                           o::Dict{Symbol,Any})
    if haskey(o, :fields)
        fields = extract_export_fields(o[:fields])
    else
        fields = default_export_fields()
    end

    check_export_fields(body, solver, fields)

    return fields
end

function get_export_fields(ms::AbstractMultibodySetup, solver::AbstractTimeSolver,
                           o::Dict{Symbol,Any})
    if haskey(o, :fields)
        fields_spec = extract_fields_spec(ms, o[:fields])
    else
        fields_spec = default_fields_spec(ms)
    end

    check_export_fields(ms, solver, fields_spec)

    return fields_spec
end

@inline function extract_export_fields(o)
    fields::Vector{Symbol} = _extract_export_fields(o)
    return fields
end

@inline _extract_export_fields(fields::Symbol) = [fields]
@inline _extract_export_fields(fields::NTuple{N,Symbol}) where {N} = [f for f in fields]
@inline _extract_export_fields(fields::Vector{Symbol}) = fields
@inline function _extract_export_fields(arbitrary_fields::T) where {T}
    fields::Vector{Symbol} = try
        convert(Vector{Symbol}, arbitrary_fields)
    catch
        msg = "wrong type of specification for keyword `fields`!\n"
        msg *= "Cannot convert input of type `$(T)` to `Vector{Symbol}`!\n"
        throw(ArgumentError(msg))
    end
    return fields
end

@inline function extract_fields_spec(ms::AbstractMultibodySetup, o)
    fields_spec::Dict{Symbol,Vector{Symbol}} = _extract_fields_spec(ms, o)
    for body_name in each_body_name(ms)
        if !haskey(fields_spec, body_name)
            fields_spec[body_name] = default_export_fields()
        end
    end
    return fields_spec
end

function _extract_fields_spec(ms::AbstractMultibodySetup, o::Dict{Symbol,T}) where {T}
    fields_spec = Dict{Symbol,Vector{Symbol}}()
    for (body_name, fields) in o
        check_if_bodyname_is_defined(ms, body_name)
        fields_spec[body_name] = extract_export_fields(fields)
    end
    return fields_spec
end

function _extract_fields_spec(ms::AbstractMultibodySetup, o)
    fields = extract_export_fields(o)
    fields_spec = Dict{Symbol,Vector{Symbol}}()
    for body_name in each_body_name(ms)
        fields_spec[body_name] = fields
    end
    return fields_spec
end

function check_export_fields(body::AbstractBody, solver::AbstractTimeSolver,
                             fields::Vector{Symbol})
    S = storage_type(body, solver)
    check_export_fields(S, fields)
    return nothing
end

function check_export_fields(ms::AbstractMultibodySetup, solver::AbstractTimeSolver,
                             fields_spec::Dict{Symbol,Vector{Symbol}})
    for body_name in each_body_name(ms)
        body = get_body(ms, body_name)
        S = storage_type(body, solver)
        fields = fields_spec[body_name]
        check_export_fields(S, fields)
    end
    return nothing
end

function check_export_fields(::Type{S}, fields::Vector{Symbol}) where {S}
    allowed_fields = point_data_fields(S)
    for f in fields
        if !in(f, allowed_fields)
            msg = "unknown point data field `:$(f)` specified for export!\n"
            msg *= "All point data fields of $S:\n"
            for allowed_name in allowed_fields
                msg *= "  - $allowed_name\n"
            end
            throw(ArgumentError(msg))
        end
    end
    return nothing
end

function get_vtk_filebase(body::AbstractBody, root::AbstractString)
    body_name = replace(string(get_name(body)), " " => "_")
    filebase = isempty(body_name) ? "timestep" : body_name * "_timestep"
    vtk_filebase::String = joinpath(root, "vtk", filebase)
    return vtk_filebase
end

function get_vtk_filebase(ms::AbstractMultibodySetup, root::AbstractString)
    vtk_filebase = Dict{Symbol,String}()
    for body_name in each_body_name(ms)
        vtk_filebase[body_name] = get_vtk_filebase(get_body(ms, body_name), root)
    end
    return vtk_filebase
end

function _export_results(options::AbstractJobOptions, chunk::AbstractBodyChunk,
                         chunk_id::Int, n_chunks::Int, n::Int, t::Float64)
    filename = get_filename(options, chunk.body_name, n)
    position = get_loc_position(chunk)
    pvtk_grid(filename, position, chunk.cells; part=chunk_id, nparts=n_chunks) do vtk
        export_fields!(vtk, chunk, options.fields)
        vtk["time", VTKFieldData()] = t
    end
    return nothing
end

@inline function get_filename(options::AbstractJobOptions, body_name::Symbol, n::Int)
    return _get_filename(options.vtk_filebase, body_name, n)
end

@inline function _get_filename(vtk_filebase::String, ::Symbol, n::Int)
    return @sprintf("%s_%06d", vtk_filebase, n)
end

@inline function _get_filename(vtk_filebase::Dict{Symbol,String}, body_name::Symbol, n::Int)
    return @sprintf("%s_%06d", vtk_filebase[body_name], n)
end

function export_fields!(vtk, chunk, fields::Vector{Symbol})
    for field in fields
        point_data = get_loc_point_data(chunk.storage, chunk.ch, field)
        vtk[string(field), VTKPointData()] = point_data
    end
    return nothing
end

function export_fields!(vtk, chunk, fields_spec::Dict{Symbol,Vector{Symbol}})
    fields = fields_spec[chunk.body_name]
    export_fields!(vtk, chunk, fields)
    return nothing
end

@inline function get_loc_position(chunk::AbstractBodyChunk)
    return @views chunk.storage.position[:, 1:chunk.ch.n_loc_points]
end
