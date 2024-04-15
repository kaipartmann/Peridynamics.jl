function storage_type(mat::AbstractMaterial, ts::AbstractTimeSolver)
    throw(MethodError(storage_type, mat, ts))
end

function init_storage(::AbstractBody{M}, ts::T, args...) where {M,T}
    msg = "storage for material $M and time solver $T not specified!\n"
    return error(msg)
end

macro storage(material, timesolver, storage)
    macrocheck_input_material(material)
    macrocheck_input_timesolver(timesolver)
    macrocheck_input_storage(storage)
    local _checks = quote
        Peridynamics.typecheck_material($(esc(material)))
        Peridynamics.typecheck_storage($(esc(storage)), $(esc(timesolver)))
    end
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
    return Expr(:block, _checks, _storage_type, _init_storage)
end

halo_read_fields(::AbstractStorage) = ()
halo_write_fields(::AbstractStorage) = ()
is_halo_field(::AbstractStorage, ::Val{F}) where {F} = false

macro halo_read_fields(storage, fields...)
    macrocheck_input_storage(storage)
    macrocheck_input_fields(fields...)
    local _checks = quote
        if !($(esc(storage)) <: Peridynamics.AbstractStorage)
            msg = Base.string($(esc(storage)), " is not a valid storage type!\n")
            throw(ArgumentError(msg))
        end
        Peridynamics.typecheck_storage_has_fields($(esc(storage)), $(Expr(:tuple, fields...)))
    end
    local _halo_read_fields = quote
        function Peridynamics.halo_read_fields(::$(esc(storage)))
            return $(Expr(:tuple, fields...))
        end
    end
    local _is_halo_fields = [
        quote
            function Peridynamics.is_halo_field(::$(esc(storage)), ::Val{$(esc(_field))})
                return true
            end
        end
        for _field in fields
    ]
    return Expr(:block, _checks, _halo_read_fields, _is_halo_fields...)
end

macro halo_write_fields(storage, fields...)
    macrocheck_input_storage(storage)
    macrocheck_input_fields(fields...)
    local _checks = quote
        if !($(esc(storage)) <: Peridynamics.AbstractStorage)
            msg = Base.string($(esc(storage)), " is not a valid storage type!\n")
            throw(ArgumentError(msg))
        end
        Peridynamics.typecheck_storage_has_fields($(esc(storage)), $(Expr(:tuple, fields...)))
    end
    local _halo_read_fields = quote
        function Peridynamics.halo_write_fields(::$(esc(storage)))
            return $(Expr(:tuple, fields...))
        end
    end
    local _is_halo_fields = [
        quote
            function Peridynamics.is_halo_field(::$(esc(storage)), ::Val{$(esc(_field))})
                return true
            end
        end
        for _field in fields
    ]
    return Expr(:block, _checks, _halo_read_fields, _is_halo_fields...)
end

function macrocheck_input_storage(storage)
    storage isa Symbol && return nothing
    (storage isa Expr && storage.head === :.) && return nothing
    msg = "argument `$storage` is not a valid storage input!\n"
    return throw(ArgumentError(msg))
end

function macrocheck_input_timesolver(timesolver)
    timesolver isa Symbol && return nothing
    (timesolver isa Expr && timesolver.head === :.) && return nothing
    msg = "argument `$timesolver` is not a valid time solver input!\n"
    return throw(ArgumentError(msg))
end

function macrocheck_input_fields(fields...)
    for field in fields
        if !(field isa QuoteNode && field.value isa Symbol)
            msg = "argument $field is not a valid field input!\n"
            throw(ArgumentError(msg))
        end
    end
    return nothing
end

function typecheck_storage(::Type{Storage}, ::Type{TimeSolver}) where {Storage,TimeSolver}
    if !(Storage <: AbstractStorage)
        msg = "$Storage is not a valid storage type!\n"
        throw(ArgumentError(msg))
    end
    # TODO: add the material type! Then users can add requirements for own materials
    req_storage_fields_timesolver(Storage, TimeSolver)
    req_storage_fields_fracture(Storage)
    return nothing
end

function typecheck_storage_has_fields(::Type{S}, fields::NTuple{N,Symbol}) where {S,N}
    storage_fields = fieldnames(S)
    for field in fields
        if !in(field, storage_fields)
            msg = "storage $(S) has no field `$(field)`!\n"
            throw(ArgumentError(msg))
        end
    end
    return nothing
end

@inline function get_point_data(s::AbstractStorage, field::Symbol)
    return point_data_field(s, Val(field))
end

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

# TODO: use required fields specifications!
function point_data_fields(::Type{S}) where {S<:AbstractStorage}
    fields = (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
              :b_ext, :damage, :n_active_bonds)
    return fields
end

function get_halo_read_fields(s::AbstractStorage)
    return Tuple(get_point_data(s, field) for field in halo_read_fields(s))
end

function get_halo_write_fields(s::AbstractStorage)
    return Tuple(get_point_data(s, field) for field in halo_write_fields(s))
end

function get_loc_point_data(s::AbstractStorage, ch::AbstractChunkHandler, field::Symbol)
    val_field = Val(field)
    is_halo_field(s, val_field) || return point_data_field(s, val_field)
    return get_loc_view(point_data_field(s, val_field), ch)
end
