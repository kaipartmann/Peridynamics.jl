function storage_type(mat::AbstractMaterial)
    return throw(InterfaceError(mat, "storage_type"))
end

function get_storage(material, solver, system)
    return throw(InterfaceError(material, "get_storage"))
end

function init_field(material, solver, system, field)
    field_value_system = init_field_system(system, field)
    isnothing(field_value_system) || return field_value_system
    field_value_solver = init_field_solver(solver, system, field)
    isnothing(field_value_solver) || return field_value_solver
    method_msg = "init_field(::$(typeof(material)), "
    method_msg *= "::$(typeof(solver)), "
    method_msg *= "::$(typeof(system)), "
    method_msg *= "::$(typeof(field)))"
    return throw(InterfaceError(M, method_msg))
end

macro storage(material, storage)
    macrocheck_input_material(material)
    macrocheck_input_storage_type(storage)
    local _checks = quote
        Peridynamics.typecheck_material($(esc(material)))
        Peridynamics.typecheck_storage($(esc(material)), $(esc(storage)))
    end
    local _storage_type = quote
        function Peridynamics.storage_type(::$(esc(material)))
            return $(esc(storage))
        end
    end
    local _constructor = quote
        function $(esc(storage))(mat::$(esc(material)),
                                solver::Peridynamics.AbstractTimeSolver,
                                system::Peridynamics.AbstractSystem)
            fields = Base.fieldnames($(esc(storage)))
            args = (Peridynamics.init_field(mat, solver, system, Val(field))
                    for field in fields)
            return $(esc(storage))(args...)
        end
    end
    local _get_storage = quote
        function Peridynamics.get_storage(mat::$(esc(material)),
                                          solver::Peridynamics.AbstractTimeSolver,
                                          system::Peridynamics.AbstractSystem)
            return $(esc(storage))(mat, solver, system)
        end
    end
    # the checks have to come last; otherwise errors cannot be tested with @test_throws ...
    return Expr(:block, _storage_type, _constructor, _get_storage, _checks)
end

# macro storagedef(material, storage)
#     macrocheck_input_material(material)
#     macrocheck_input_storage_struct(storage)

#     # extract the fields of the struct, annotated
#     local _point_fields = Vector{Expr}()
#     local _htl_fields = Vector{Expr}()
#     local _lth_fields = Vector{Expr}()
#     local _constructor = quote
#         function $(esc(storage))(mat::$(esc(material)),
#                                  solver::Peridynamics.AbstractTimeSolver,
#                                  system::Peridynamics.AbstractSystem)
#             fields = Base.fieldnames($(esc(storage)))
#             args = (Peridynamics.init_field(mat, solver, system, Val(field))
#                     for field in fields)
#             return $(esc(storage))(args...)
#         end
#     end
#     local _storage_type = quote
#         function Peridynamics.storage_type(::$(esc(material)))
#             return $(esc(storage))
#         end
#     end
#     local _get_storage = quote
#         function Peridynamics.get_storage(mat::$(esc(material)),
#                                           solver::Peridynamics.AbstractTimeSolver,
#                                           system::Peridynamics.AbstractSystem)
#             return $(esc(storage))(mat, solver, system)
#         end
#     end
#     local _checks = quote
#         Peridynamics.typecheck_material($(esc(material)))
#         Peridynamics.typecheck_storage($(esc(material)), $(esc(storage)))
#     end
#     local _exprs = Expr(:block, _storage_type, _constructor, _get_storage, _checks)
#     return _exprs
# end

# macro pointfield(expr)
#     return nothing
# end

# macro htlfield(expr)
#     return nothing
# end

# macro lthfield(expr)
#     return nothing
# end

loc_to_halo_fields(::AbstractStorage) = ()
halo_to_loc_fields(::AbstractStorage) = ()
is_halo_field(::AbstractStorage, ::Val{F}) where {F} = false

macro loc_to_halo_fields(storage, fields...)
    macrocheck_input_storage_type(storage)
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
    # the checks have to come last; otherwise errors cannot be tested with @test_throws ...
    return Expr(:block, _loc_to_halo_fields, _is_halo_fields..., _checks)
end

macro halo_to_loc_fields(storage, fields...)
    macrocheck_input_storage_type(storage)
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
    # the checks have to come last; otherwise errors cannot be tested with @test_throws ...
    return Expr(:block, _loc_to_halo_fields, _is_halo_fields..., _checks)
end

macro halo_fields(storage, fields...)
    macrocheck_input_storage_type(storage)
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
    # the checks have to come last; otherwise errors cannot be tested with @test_throws ...
    return Expr(:block, _loc_to_halo_fields, _is_halo_fields..., _checks)
end

function macrocheck_input_storage_type(storage)
    storage isa Symbol && return nothing
    (storage isa Expr && storage.head === :.) && return nothing
    msg = "argument `$storage` is not a valid storage input!\n"
    return throw(ArgumentError(msg))
end

function macrocheck_input_storage_struct(storage)
    if !(storage isa Expr && storage.head === :struct)
        msg = "specified input is not a valid point parameter struct expression!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function macrocheck_input_timesolver(timesolver)
    timesolver isa Symbol && return nothing
    (timesolver isa Expr && timesolver.head === :.) && return nothing
    msg = "argument `$timesolver` is not a valid time solver input!\n"
    return throw(ArgumentError(msg))
end

function macrocheck_input_fields(fields...)
    for field in fields
        macrocheck_input_field(field)
    end
    return nothing
end

function macrocheck_input_field(field)
    if !(field isa QuoteNode && field.value isa Symbol)
        msg = "argument $field is not a valid field input!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function typecheck_storage(::Type{Material}, ::Type{Storage}) where {Material,Storage}
    typecheck_is_storage(Storage)
    typecheck_req_fields_missing(Storage, required_fields(Material))
    return nothing
end

function typecheck_storage(material, storage)
    return throw(ArgumentError("$storage is not a valid storage type!\n"))
end

function typecheck_is_storage(::Type{Storage}) where {Storage}
    if !(Storage <: AbstractStorage)
        msg = "$Storage is not a valid storage type!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function typecheck_is_storage(storage)
    return throw(ArgumentError("$storage is not a valid storage type!\n"))
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

function required_fields(::Type{Material}) where {Material<:AbstractMaterial}
    fields = (required_fields_timesolvers()..., required_fields_fracture(Material)...)
    return fields
end

@inline function get_point_data(s::AbstractStorage, field::Symbol)
    return point_data_field(s, Val(field))
end

function point_data_field(::S, ::Val{Field}) where {S,Field}
    return throw(InterfaceError(S, "point_data_field(::$S, ::Val{$Field})"))
end

# macro point_data_field(storage, field)
#     macrocheck_input_storage_type(storage)
#     macrocheck_input_field(field)
#     local _field_symb = field isa QuoteNode ? field : QuoteNode(field)
#     local _pdf_func = quote
#         function Peridynamics.point_data_field(storage::$(esc(storage)),
#                                                ::Base.Val{$(_field_symb)})
#             return Base.getfield(storage, $(_field_symb))
#         end
#     end
#     return _pdf_func
# end

point_data_field(s::AbstractStorage, ::Val{:position}) = getfield(s, :position)
point_data_field(s::AbstractStorage, ::Val{:displacement}) = getfield(s, :displacement)
point_data_field(s::AbstractStorage, ::Val{:velocity}) = getfield(s, :velocity)
point_data_field(s::AbstractStorage, ::Val{:velocity_half}) = getfield(s, :velocity_half)
point_data_field(s::AbstractStorage, ::Val{:acceleration}) = getfield(s, :acceleration)
point_data_field(s::AbstractStorage, ::Val{:b_int}) = getfield(s, :b_int)
point_data_field(s::AbstractStorage, ::Val{:b_ext}) = getfield(s, :b_ext)
point_data_field(s::AbstractStorage, ::Val{:damage}) = getfield(s, :damage)
point_data_field(s::AbstractStorage, ::Val{:n_active_bonds}) = getfield(s, :n_active_bonds)
point_data_field(s::AbstractStorage, ::Val{:n_active_one_nis}) = getfield(s, :n_active_one_nis)

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
