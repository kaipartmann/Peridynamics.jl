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

    S = storage_type(body, solver)
    check_export_fields(S, fields)

    return fields
end

function get_export_fields(ms::AbstractMultibodySetup, solver::AbstractTimeSolver,
                           o::Dict{Symbol,Any})
    if haskey(o, :fields)
        fields_spec = extract_fields_spec(ms, o[:fields])
    else
        fields_spec = default_fields_spec(ms)
    end

    return fields_spec
end

@inline extract_export_fields(fields::Symbol) = [fields]
@inline extract_export_fields(fields::NTuple{N,Symbol}) where {N} = [f for f in fields]
@inline extract_export_fields(fields::Vector{Symbol}) = fields
@inline function extract_export_fields(arbitrary_fields::T) where {T}
    fields::Vector{Symbol} = try
        convert(Vector{Symbol}, arbitrary_fields)
    catch
        msg = "wrong type of specification for keyword `fields`!\n"
        msg *= "Cannot convert input of type `$(T)` to `Vector{Symbol}`!\n"
        throw(ArgumentError(msg))
    end
    return fields
end

function extract_fields_spec(ms::AbstractMultibodySetup, o::Dict{Symbol,T}) where {T}
    fields_spec = Dict{Symbol,Vector{Symbol}}()
    for (body_name, fields) in o
        check_if_bodyname_is_defined(ms, body_name)
        fields_spec[body_name] = extract_export_fields(fields)
    end
    return fields_spec
end

function extract_fields_spec(ms::AbstractMultibodySetup, o)
    return extract_export_fields(o)
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

function _export_results(options::AbstractJobOptions, b::AbstractBodyChunk, chunk_id::Int,
                         n_chunks::Int, prefix::AbstractString, n::Int, t::Float64)
    filename = joinpath(options.vtk, @sprintf("%s_timestep_%06d", prefix, n))
    position = get_loc_position(b)
    pvtk_grid(filename, position, b.cells; part=chunk_id, nparts=n_chunks) do vtk
        export_fields!(vtk, b, options.fields)
        vtk["time", VTKFieldData()] = t
    end
    return nothing
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

@inline function get_loc_position(b::AbstractBodyChunk)
    return @views b.storage.position[:, 1:b.ch.n_loc_points]
end
