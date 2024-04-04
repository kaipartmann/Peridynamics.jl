
function get_export_fields(::Type{S}, o::Dict{Symbol,Any}) where {S}
    local _fields::NTuple{N,Symbol} where {N}
    if haskey(o, :fields)
        _fields = o[:fields]
    else
        _fields = DEFAULT_EXPORT_FIELDS
    end

    fields = [field for field in _fields]
    check_export_fields(S, fields)

    return fields
end

function check_export_fields(::Type{S}, fields::Vector{Symbol}) where {S}
    allowed_fields = fieldnames(S)
    for f in fields
        if !in(f, allowed_fields)
            msg = "unknown field $(name) specified for export!\n"
            msg *= "Allowed fields for $S:\n"
            for allowed_name in fieldnames(S)
                msg *= "  - $allowed_name\n"
            end
            throw(ArgumentError(msg))
        end
    end
    return nothing
end


function _export_results(b::AbstractBodyChunk, chunk_id::Int, n_chunks::Int,
                         options::AbstractOptions, n::Int, t::Float64)
    filename = joinpath(options.vtk, @sprintf("timestep_%06d", n))
    position = get_loc_position(b)
    pvtk_grid(filename, position, b.cells; part=chunk_id, nparts=n_chunks) do vtk
        for field in options.fields
            vtk[string(field), VTKPointData()] = get_storage_field(b.storage, field)
        end
        vtk["time", VTKFieldData()] = t
    end
    return nothing
end

@inline function get_loc_position(b::AbstractBodyChunk)
    return @views b.storage.position[:, 1:b.ch.n_loc_points]
end
