@inline function get_storage_field(s::AbstractStorage, field::Symbol)
    return storage_field(s, Val(field))
end

storage_field(s::AbstractStorage, ::Val{:position}) = getfield(s, :position)
storage_field(s::AbstractStorage, ::Val{:displacement}) = getfield(s, :displacement)
storage_field(s::AbstractStorage, ::Val{:velocity}) = getfield(s, :velocity)
storage_field(s::AbstractStorage, ::Val{:velocity_half}) = getfield(s, :velocity_half)
storage_field(s::AbstractStorage, ::Val{:acceleration}) = getfield(s, :acceleration)
storage_field(s::AbstractStorage, ::Val{:b_int}) = getfield(s, :b_int)
storage_field(s::AbstractStorage, ::Val{:b_ext}) = getfield(s, :b_ext)
storage_field(s::AbstractStorage, ::Val{:damage}) = getfield(s, :damage)

function storage_field(s::AbstractStorage, ::Val{field}) where {field}
    return getfield(s, field)
end

function get_halo_read_fields(s::AbstractStorage)
    return Tuple(get_storage_field(s, field) for field in halo_read_fields(s))
end

function get_halo_write_fields(s::AbstractStorage)
    return Tuple(get_storage_field(s, field) for field in halo_write_fields(s))
end
