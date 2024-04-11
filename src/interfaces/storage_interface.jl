
@inline function get_point_data(s::AbstractStorage, field::Symbol)
    return point_data_field(s, Val(field))
end

# const DEFAULT_POINT_DATA_FIELDS = (:position, :displacement, :velocity, :velocity_half,
#                                    :acceleration, :b_int, :b_ext, :damage, :n_active_bonds)

# TODO: create these functions with a macro
point_data_field(s::AbstractStorage, ::Val{:position}) = getfield(s, :position)
point_data_field(s::AbstractStorage, ::Val{:displacement}) = getfield(s, :displacement)
point_data_field(s::AbstractStorage, ::Val{:velocity}) = getfield(s, :velocity)
point_data_field(s::AbstractStorage, ::Val{:velocity_half}) = getfield(s, :velocity_half)
point_data_field(s::AbstractStorage, ::Val{:acceleration}) = getfield(s, :acceleration)
point_data_field(s::AbstractStorage, ::Val{:b_int}) = getfield(s, :b_int)
point_data_field(s::AbstractStorage, ::Val{:b_ext}) = getfield(s, :b_ext)
point_data_field(s::AbstractStorage, ::Val{:damage}) = getfield(s, :damage)
point_data_field(s::AbstractStorage, ::Val{:n_active_bonds}) = getfield(s, :n_active_bonds)

function point_data_field(::AbstractStorage, ::Val{field}) where {field}
    return throw(ArgumentError("field $field is not a valid point data field!\n"))
end

function get_halo_read_fields(s::AbstractStorage)
    return Tuple(get_point_data(s, field) for field in halo_read_fields(s))
end

function get_halo_write_fields(s::AbstractStorage)
    return Tuple(get_point_data(s, field) for field in halo_write_fields(s))
end

function point_data_fields(s::Type{S}) where {S<:AbstractStorage}
    return throw(MethodError(point_data_fields, s))
end

function get_loc_point_data(s::AbstractStorage, ch::AbstractChunkHandler, field::Symbol)
    val_field = Val(field)
    is_halo_field(s, val_field) || return point_data_field(s, val_field)
    return get_loc_view(point_data_field(s, val_field), ch)
end

##----
function init_storage(::AbstractBody{M}, ts::T, args...) where {M,T}
    msg = "storage for material $M and time solver $T not specified!\n"
    return error(msg)
end

macro storage(material, timesolver, storage)
    local _storage_type = quote
        function Peridynamics.storage_type(::$(esc(material)), ::$(esc(timesolver)))
            return $(esc(storage))
        end
    end
    local _init_storage = quote
        function Peridynamics.init_storage(body::Peridynamics.AbstractBody{$(esc(material))},
                                           ts::$(esc(timesolver)), args...)
            return $(esc(storage))(body, ts, args...)
        end
    end
    return Expr(:block, _storage_type, _init_storage)
end
