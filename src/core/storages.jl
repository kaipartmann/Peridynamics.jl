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

macro storagedef(material, storage)
    macrocheck_input_material(material)
    macrocheck_input_storage_struct(storage)

    local _expr_extraction = get_storage_structdef(storage)
    local _storage_struct, _point_fields, _htl_fields, _lth_fields = _expr_extraction
    # local _point_field_names = Tuple(f.args[1] for f in _point_fields)
    local _htl_field_names = Tuple(QuoteNode(f.args[1]) for f in _htl_fields)
    local _lth_field_names = Tuple(QuoteNode(f.args[1]) for f in _lth_fields)
    local _halo_field_names = (_htl_field_names..., _lth_field_names...)
    local _storage_type = extract_storage_type(_storage_struct)
    local _constructor = quote
        function $(esc(_storage_type))(mat::$(esc(material)),
                                 solver::Peridynamics.AbstractTimeSolver,
                                 system::Peridynamics.AbstractSystem)
            fields = Base.fieldnames($(esc(_storage_type)))
            args = (Peridynamics.init_field(mat, solver, system, Val(field))
                    for field in fields)
            return $(esc(_storage_type))(args...)
        end
    end
    local _storage_type_function = quote
        function Peridynamics.storage_type(::$(esc(material)))
            return $(esc(_storage_type))
        end
    end
    local _get_storage = quote
        function Peridynamics.get_storage(mat::$(esc(material)),
                                          solver::Peridynamics.AbstractTimeSolver,
                                          system::Peridynamics.AbstractSystem)
            return $(esc(_storage_type))(mat, solver, system)
        end
    end
    local _halo_to_loc_fields = quote
        function Peridynamics.halo_to_local_fields(::$(_storage_type))
            return $(Expr(:tuple, _htl_field_names...))
        end
    end
    local _loc_to_halo_fields = quote
        function Peridynamics.loc_to_halo_fields(::$(_storage_type))
            return $(Expr(:tuple, _lth_field_names...))
        end
    end
    local _is_halo_fields = [quote
                                 function Peridynamics.is_halo_field(::$(_storage_type),
                                                                     ::Val{$(_field)})
                                     return true
                                 end
                             end
                             for _field in _halo_field_names]
    local _checks = quote
        Peridynamics.typecheck_material($(esc(material)))
        Peridynamics.typecheck_storage($(esc(material)), $(esc(_storage_type)))
    end
    local _exprs = Expr(:block, _storage_struct, _constructor, _storage_type_function,
                        _get_storage, _halo_to_loc_fields, _loc_to_halo_fields,
                        _is_halo_fields, _checks)
    return _exprs
end

function get_storage_structdef(storage_expr)
    # extract the fields of the struct
    fields_exprs = Vector{Expr}()

    # indices of special annotated fields in `_fields`
    point_fields_idxs = Vector{Int}()
    htl_fields_idxs = Vector{Int}()
    lth_fields_idxs = Vector{Int}()

    field_idx = 0
    for field in storage_expr.args[3].args
        field isa Expr || continue
        # just a normal data field (bond data or arbitrary data)
        if field.head === :(::)
            field_idx += 1
            push!(fields_exprs, field)
        # a point data field
        elseif field.head === :macrocall && field.args[1] === Symbol("@pointfield")
            annotated_field = field.args[3]
            if annotated_field isa Expr && annotated_field.head == :(::)
                field_idx += 1
                push!(fields_exprs, annotated_field)
                push!(point_fields_idxs, field_idx)
            else
                msg = "unexpected or untyped field: $annotated_field\n"
                throw(ArgumentError(msg))
            end
        elseif field.head === :macrocall && field.args[1] === Symbol("@htlfield")
            annotated_field = field.args[3]
            if annotated_field isa Expr && annotated_field.head == :(::)
                field_idx += 1
                push!(fields_exprs, annotated_field)
                push!(htl_fields_idxs, field_idx)
            else
                msg = "unexpected or untyped field: $annotated_field\n"
                throw(ArgumentError(msg))
            end
        elseif field.head === :macrocall && field.args[1] === Symbol("@lthfield")
            annotated_field = field.args[3]
            if annotated_field isa Expr && annotated_field.head == :(::)
                field_idx += 1
                push!(fields_exprs, annotated_field)
                push!(lth_fields_idxs, field_idx)
            else
                msg = "unexpected or untyped field: $annotated_field\n"
                throw(ArgumentError(msg))
            end
        else
            msg = "unexpected or untyped field: $field\n"
            throw(ArgumentError(msg))
        end
    end
    fields_block = Expr(:block, fields_exprs...)
    storage_struct = Expr(:struct, storage_expr.args[1], storage_expr.args[2], fields_block)
    point_fields_exprs = [fields_exprs[i] for i in point_fields_idxs]
    htl_fields_exprs = [fields_exprs[i] for i in htl_fields_idxs]
    lth_fields_exprs = [fields_exprs[i] for i in lth_fields_idxs]
    return storage_struct, point_fields_exprs, htl_fields_exprs, lth_fields_exprs
end

function extract_storage_type(storage_expr)
    header = storage_expr.args[2]
    if header.head === :(<:) && header.args[1] isa Symbol
        return header.args[1]
    else
        error("not yet supported storage header: `$header`\n")
    end
end

macro pointfield(expr)
    return nothing
end

macro htlfield(expr)
    return nothing
end

macro lthfield(expr)
    return nothing
end

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
