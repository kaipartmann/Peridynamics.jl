function storage_type(mat::AbstractMaterial)
    return throw(InterfaceError(mat, "storage_type"))
end

function get_storage(material, solver, system)
    return throw(InterfaceError(material, "get_storage"))
end

function init_field(material, solver, system, field)
    method = "init_field(::$(typeof(material)), "
    method *= "::$(typeof(solver)), "
    method *= "::$(typeof(system)), "
    method *= "::$(typeof(field)))"
    return throw(InterfaceError(M, method))
end

macro storage(material, storage)
    macrocheck_input_material(material)
    macrocheck_input_storage(storage)
    local _checks = quote
        Peridynamics.typecheck_material($(esc(material)))
        Peridynamics.typecheck_storage($(esc(storage)))
    end
    local _storage_type = quote
        function Peridynamics.storage_type(::$(esc(material)))
            return $(esc(storage))
        end
    end
    local _constructor = quote
        function $(esc(storage))(mat::$(esc(material)), solver::AbstractTimeSolver,
                                 system::AbstractSystem)
            fields = fieldnames($(esc(storage)))
            args = (init_field(mat, solver, system, Val(field)) for field in fields)
            return $(esc(storage))(args...)
        end
    end
    local _get_storage = quote
        function Peridynamics.get_storage(mat::$(esc(material)), solver::AbstractTimeSolver,
                                          system::AbstractSystem)
            return $(esc(storage))(mat, solver, system)
        end
    end
    return Expr(:block, _checks, _storage_type, _constructor, _get_storage)
end

loc_to_halo_fields(::AbstractStorage) = ()
halo_to_loc_fields(::AbstractStorage) = ()
is_halo_field(::AbstractStorage, ::Val{F}) where {F} = false

macro loc_to_halo_fields(storage, fields...)
    macrocheck_input_storage(storage)
    macrocheck_input_fields(fields...)
    local _checks = quote
        Peridynamics.typecheck_is_storage($(esc(storage)))
        Peridynamics.typecheck_storage_fields($(esc(storage)), $(Expr(:tuple, fields...)))
    end
    local _loc_to_halo_fields = quote
        function Peridynamics.loc_to_halo_fields(::$(esc(storage)))
            return $(Expr(:tuple, fields...))
        end
    end
    local _is_halo_fields = [quote
                                 function Peridynamics.is_halo_field(::$(esc(storage)),
                                                                     ::Val{$(esc(_field))})
                                     return true
                                 end
                             end
                             for _field in fields]
    return Expr(:block, _checks, _loc_to_halo_fields, _is_halo_fields...)
end

macro halo_to_loc_fields(storage, fields...)
    macrocheck_input_storage(storage)
    macrocheck_input_fields(fields...)
    local _checks = quote
        Peridynamics.typecheck_is_storage($(esc(storage)))
        Peridynamics.typecheck_storage_fields($(esc(storage)), $(Expr(:tuple, fields...)))
    end
    local _loc_to_halo_fields = quote
        function Peridynamics.halo_to_loc_fields(::$(esc(storage)))
            return $(Expr(:tuple, fields...))
        end
    end
    local _is_halo_fields = [quote
                                 function Peridynamics.is_halo_field(::$(esc(storage)),
                                                                     ::Val{$(esc(_field))})
                                     return true
                                 end
                             end
                             for _field in fields]
    return Expr(:block, _checks, _loc_to_halo_fields, _is_halo_fields...)
end

macro halo_fields(storage, fields...)
    macrocheck_input_storage(storage)
    macrocheck_input_fields(fields...)
    local _checks = quote
        Peridynamics.typecheck_is_storage($(esc(storage)))
        Peridynamics.typecheck_storage_fields($(esc(storage)), $(Expr(:tuple, fields...)))
    end
    local _is_halo_fields = [quote
                                 function Peridynamics.is_halo_field(::$(esc(storage)),
                                                                     ::Val{$(esc(_field))})
                                     return true
                                 end
                             end
                             for _field in fields]
    return Expr(:block, _checks, _loc_to_halo_fields, _is_halo_fields...)
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

function typecheck_is_storage(::Type{Storage}) where {Storage}
    if !(Storage <: AbstractStorage)
        msg = "$Storage is not a valid storage type!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function typecheck_storage(::Type{Storage}) where {Storage}
    typecheck_is_storage(Storage)
    typecheck_req_fields_missing(Storage, required_fields_timesolvers())
    # TODO: this check would currently be system dependent!
    # However, this is still not optimal. There has to be some additional type for different
    # damage models...
    # typecheck_req_fields_missing(Storage, required_fields_fracture())
    return nothing
end

function typecheck_storage_fields(::Type{S}, fields::NTuple{N,Symbol}) where {S,N}
    storage_fields = fieldnames(S)
    for field in fields
        if !in(field, storage_fields)
            msg = "storage $(S) has no field `$(field)`!\n"
            throw(ArgumentError(msg))
        end
    end
    return nothing
end

function typecheck_req_fields_missing(::Type{Storage}, req_fields) where {Storage}
    storage_fields = fieldnames(Storage)
    for req_field in req_fields
        if !in(req_field, storage_fields)
            error("required field $req_field not found in $(Storage)!\n")
        end
    end
    return false
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
point_data_field(s::AbstractStorage, ::Val{:n_active_one_nis}) = getfield(s, :n_active_bonds)

function point_data_field(::S, ::Val{Field}) where {S,Field}
    return throw(InterfaceError(S, "point_data_field(::$S, ::Val{$Field})"))
end

# macro point_data_field(fields)
#     if fields isa Expr && fields.head === :tuple
#         local pdfields = fields.args
#     else
#         local pdfields = eval(fields)
#     end
#     local exprs = Expr[]
#     for _field in pdfields
#         if _field isa QuoteNode
#             field = _field
#         elseif _field isa Symbol
#             field = QuoteNode(_field)
#         else
#             msg = "input $(field) is not a QuoteNode or a Symbol!\n"
#             throw(ArgumentError(msg))
#         end
#         local _get_point_data_function = quote
#             function Peridynamics.point_data_field(s::Peridynamics.AbstractStorage,
#                                                    ::Base.Val{$(field)})
#                 return Base.getfield(s, $(field))
#             end
#         end
#         push!(exprs, _get_point_data_function)
#     end
#     return Expr(:block, exprs...)
# end

# @point_data_field required_fields_timesolvers()

# TODO: use required fields specifications!
function point_data_fields(::Type{S}) where {S<:AbstractStorage}
    fields = (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
              :b_ext, :damage, :n_active_bonds)
    return fields
end

function get_loc_to_halo_fields(s::AbstractStorage)
    return Tuple(get_point_data(s, field) for field in loc_to_halo_fields(s))
end

function get_halo_to_loc_fields(s::AbstractStorage)
    return Tuple(get_point_data(s, field) for field in halo_to_loc_fields(s))
end

function get_halo_fields(s::S) where {S<:AbstractStorage}
    return Tuple(get_point_data(s, f) for f in fieldnames(S) if is_halo_field(s, Val(f)))
end

function get_halo_fieldnames(s::S) where {S<:AbstractStorage}
    return Tuple(f for f in fieldnames(S) if is_halo_field(s, Val(f)))
end

function get_n_halo_fields(s::S) where {S<:AbstractStorage}
    n_halo_fields = 0
    for f in fieldnames(S)
        if is_halo_field(s, Val(f))
            n_halo_fields += 1
        end
    end
    return n_halo_fields
end

function get_loc_point_data(storage::AbstractStorage, system::AbstractSystem, field::Symbol)
    val_field = Val(field)
    is_halo_field(storage, val_field) || return point_data_field(storage, val_field)
    return get_loc_view(point_data_field(storage, val_field), system)
end

@inline function get_coordinates(s, i)
    return SVector{3}(s.position[1, i], s.position[2, i], s.position[3, i])
end

@inline function get_coordinates_diff(s, i, j)
    return SVector{3}(s.position[1, j] - s.position[1, i],
                      s.position[2, j] - s.position[2, i],
                      s.position[3, j] - s.position[3, i])
end

@inline function get_diff(a, i, j)
    return SVector{3}(a[1, j] - a[1, i], a[2, j] - a[2, i], a[3, j] - a[3, i])
end

@inline function update_add_b_int!(storage::AbstractStorage, i::Int, b::SVector{3})
    storage.b_int[1, i] += b[1]
    storage.b_int[2, i] += b[2]
    storage.b_int[3, i] += b[3]
    return nothing
end
