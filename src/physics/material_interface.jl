
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

@inline function get_storage_field(s::AbstractStorage, field::Symbol)
    return storage_field(s, Val{field}())
end

storage_field(s::AbstractStorage, ::Val{:position}) = getfield(s, :position)
storage_field(s::AbstractStorage, ::Val{:displacement}) = getfield(s, :displacement)
storage_field(s::AbstractStorage, ::Val{:velocity}) = getfield(s, :velocity)
storage_field(s::AbstractStorage, ::Val{:velocity_half}) = getfield(s, :velocity_half)
storage_field(s::AbstractStorage, ::Val{:acceleration}) = getfield(s, :acceleration)
storage_field(s::AbstractStorage, ::Val{:b_int}) = getfield(s, :b_int)
storage_field(s::AbstractStorage, ::Val{:b_ext}) = getfield(s, :b_ext)
storage_field(s::AbstractStorage, ::Val{:damage}) = getfield(s, :damage)
function storage_field(::S, ::Val{F}) where {S,F}
    throw(ArgumentError("field $F no known in $(S)\n"))
end
