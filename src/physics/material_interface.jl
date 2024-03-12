
function point_param_type(mat::AbstractMaterial)
    throw(MethodError(point_param_type, mat))
end

function allowed_material_kwargs(mat::AbstractMaterial)
    throw(MethodError(allowed_material_kwargs, mat))
end

function get_point_params(mat::AbstractMaterial, ::Dict{Symbol,Any})
    throw(MethodError(get_point_params, mat))
end

function discretization_type(mat::AbstractMaterial)
    throw(MethodError(discretization_type, mat))
end

function storage_type(mat::AbstractMaterial, ts::AbstractTimeSolver)
    throw(MethodError(storage_type, mat, ts))
end

function get_halo_read_fields(s::AbstractStorage)
    throw(MethodError(get_halo_read_fields, s))
end

function get_halo_write_fields(s::AbstractStorage)
    throw(MethodError(get_halo_write_fields, s))
end

function get_storage_field(s::S, field::Symbol) where {S<:AbstractStorage}
    if field === :position
        return getfield(s, :position)
    elseif field === :displacement
        return getfield(s, :displacement)
    elseif field === :velocity
        return getfield(s, :velocity)
    elseif field === :velocity_half
        return getfield(s, :velocity_half)
    elseif field === :acceleration
        return getfield(s, :acceleration)
    elseif field === :b_int
        return getfield(s, :b_int)
    elseif field === :b_ext
        return getfield(s, :b_ext)
    elseif field === :damage
        return getfield(s, :damage)
    else
        return get_custom_field(s, field)
    end
end

function get_custom_field(::S, field::Symbol) where {S<:AbstractStorage}
    return throw(ArgumentError("field $field no known in $(S)\n"))
end
