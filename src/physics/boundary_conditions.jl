"""
    SingleDimBC{F}

$(internal_api_warning())

Type for a boundary condition in a single dimension for a peridynamic simulation.

# Type Parameters

- `F<:Function`: Time dependent function which describes the boundary condition.

# Fields

- `fun::F`: Time dependent function which describes the boundary condition.
- `field::Symbol`: Field of the condition (e.g. velocity, force density).
- `point_set::Symbol`: Point set on which the condition is applied.
- `dim::UInt8`: Dimension in which the condition is applied.
"""
struct SingleDimBC{F<:Function} <: AbstractCondition
    fun::F
    field::Symbol
    point_set::Symbol
    dim::UInt8
end

function Base.show(io::IO, @nospecialize(bc::SingleDimBC))
    print(io, "BC on ", field_to_name(bc.field), ": ")
    print(io, msg_fields_inline(bc, (:point_set, :dim)))
    return nothing
end

@inline function (b::SingleDimBC{F})(t::Float64) where {F}
    value::Float64 = b.fun(t)
    return value
end

function apply_bc!(s::AbstractStorage, psets::Dict{Symbol,Vector{Int}},
                   bc::SingleDimBC{F}, time::Float64) where {F}
    value = bc(time)
    isfinite(value) || return nothing
    apply_sdbc!(get_point_data(s, bc.field), value, bc.dim, psets[bc.point_set])
    return nothing
end

@inline function apply_sdbc!(field::Matrix{Float64}, value::Float64, dim::UInt8,
                             point_ids::Vector{Int})
    @simd for i in point_ids
        @inbounds field[dim, i] = value
    end
    return nothing
end

"""
    PosDepSingleDimBC{F}

$(internal_api_warning())

Type for a position dependent boundary condition in a single dimension for a peridynamic
simulation.

# Type Parameters

- `F<:Function`: Position and time dependent function which describes the boundary
    condition.

# Fields

- `fun::F`: Position and time dependent function which describes the boundary condition.
- `field::Symbol`: Field of the condition (e.g. velocity, force density).
- `point_set::Symbol`: Point set on which the condition is applied.
- `dim::UInt8`: Dimension in which the condition is applied.
"""
struct PosDepSingleDimBC{F<:Function} <: AbstractCondition
    fun::F
    field::Symbol
    point_set::Symbol
    dim::UInt8
end

function Base.show(io::IO, @nospecialize(bc::PosDepSingleDimBC))
    print(io, "Pos.-dep. BC on ", field_to_name(bc.field), ": ")
    print(io, msg_fields_inline(bc, (:point_set, :dim)))
    return nothing
end

@inline function (b::PosDepSingleDimBC{F})(p::AbstractVector, t::Float64) where {F}
    value::Float64 = b.fun(p, t)
    return value
end

function apply_bc!(s::AbstractStorage, psets::Dict{Symbol,Vector{Int}},
                   bc::PosDepSingleDimBC{F}, position::Matrix{Float64},
                   time::Float64) where {F}
    apply_pdsdbc!(get_point_data(s, bc.field), position, bc, psets[bc.point_set], time)
    return nothing
end

@inline function apply_pdsdbc!(field::Matrix{Float64}, pos::Matrix{Float64},
                               bc::PosDepSingleDimBC{F}, point_ids::Vector{Int},
                               t::Float64) where {F}
    @simd for i in point_ids
        value = bc(SVector{3}(pos[1, i], pos[2, i], pos[3, i]), t)
        if isfinite(value)
            field[bc.dim, i] = value
        end
    end
    return nothing
end

"""
    PosSingleDimBC{F}

$(internal_api_warning())

Type for a position dependent boundary condition in a single dimension for a peridynamic
simulation.

# Type Parameters

- `F<:Function`: A position dependent function which describes the boundary
    condition, not time dependent.

# Fields

- `fun::F`: Position dependent function which describes the boundary condition.
- `field::Symbol`: Field of the condition (e.g. velocity, force density).
- `point_set::Symbol`: Point set on which the condition is applied.
- `dim::UInt8`: Dimension in which the condition is applied.
"""
struct PosSingleDimBC{F<:Function} <: AbstractCondition
    fun::F
    field::Symbol
    point_set::Symbol
    dim::UInt8
end

function Base.show(io::IO, @nospecialize(bc::PosSingleDimBC))
    print(io, "Pos. BC on ", field_to_name(bc.field), ": ")
    print(io, msg_fields_inline(bc, (:point_set, :dim)))
    return nothing
end

@inline function (b::PosSingleDimBC{F})(p::AbstractVector) where {F}
    value::Float64 = b.fun(p)
    return value
end

function apply_bc!(s::AbstractStorage, psets::Dict{Symbol,Vector{Int}},
                   bc::PosSingleDimBC{F}, position::Matrix{Float64}, β::Float64) where {F}
    apply_psdbc!(get_point_data(s, bc.field), position, bc, psets[bc.point_set], β)
    return nothing
end

@inline function apply_psdbc!(field::Matrix{Float64}, pos::Matrix{Float64},
                              bc::PosSingleDimBC{F}, point_ids::Vector{Int},
                              β::Float64) where {F}
    @simd for i in point_ids
        value = β * bc(SVector{3}(pos[1, i], pos[2, i], pos[3, i]))
        if isfinite(value)
            field[bc.dim, i] = value
        end
    end
    return nothing
end

function constrained_dofs!(constrained_dofs::Vector{Int}, system::AbstractSystem,
                           bc::PosSingleDimBC, loc_point_sets::Dict{Symbol,Vector{Int}})
    if bc.field === :displacement
        points = loc_point_sets[bc.point_set]
        for (dof, d, _) in each_dof_idx(system, points)
            if d == bc.dim
                push!(constrained_dofs, dof)
            end
        end
    end
    return nothing
end

struct DataBC <: AbstractCondition
    data::Matrix{Float64}
    field::Symbol
    point_set::Symbol
    dims::Vector{UInt8}
end

function Base.show(io::IO, @nospecialize(bc::DataBC))
    print(io, "Data BC on ", field_to_name(bc.field), ": ")
    print(io, msg_fields_inline(bc, (:point_set, :dims)))
    return nothing
end

function apply_bc!(s::AbstractStorage, psets::Dict{Symbol,Vector{Int}}, bc::DataBC,
                   all_point_ids::Vector{Int})
    apply_databc!(get_point_data(s, bc.field), bc, psets[bc.point_set], all_point_ids)
    return nothing
end

@inline function apply_databc!(field::Matrix{Float64}, bc::DataBC, point_ids::Vector{Int},
                               all_point_ids::Vector{Int})
    for i in point_ids
        for j in eachindex(bc.dims)
            value = bc.data[j, all_point_ids[i]]
            if isfinite(value)
                dim = bc.dims[j]
                field[dim, i] = value
            end
        end
    end
    return nothing
end

function apply_boundary_conditions!(chunk::AbstractBodyChunk, t)
    (; storage, condhandler) = chunk
    (; single_dim_bcs, posdep_single_dim_bcs, data_bcs) = condhandler
    (; loc_point_sets) = condhandler
    for bc in single_dim_bcs
        apply_bc!(storage, loc_point_sets, bc, t)
    end
    for bc in posdep_single_dim_bcs
        apply_bc!(storage, loc_point_sets, bc, chunk.system.position, t)
    end
    for bc in data_bcs
        apply_bc!(storage, loc_point_sets, bc, chunk.system.chunk_handler.point_ids)
    end
    return nothing
end

# Only for the Newton-Raphson solver, incremental displacement BCs
function apply_incr_boundary_conditions!(chunk::AbstractBodyChunk, β)
    (; storage, condhandler) = chunk
    (; pos_single_dim_bcs, loc_point_sets) = condhandler
    for bc in pos_single_dim_bcs
        apply_bc!(storage, loc_point_sets, bc, chunk.system.position, β)
    end
    return nothing
end

function constrained_dofs!(constrained_dofs::Vector{Int}, body::AbstractBody,
                           system::AbstractSystem, loc_point_sets::Dict{Symbol,Vector{Int}})
    (; pos_single_dim_bcs) = body # for now, only single dim bcs constrain dofs
    for bc in pos_single_dim_bcs
        constrained_dofs!(constrained_dofs, system, bc, loc_point_sets)
    end
    return nothing
end

function wrong_dim_err_msg(dim)
    msg = "unknown dimension `$(dim)`!\n"
    msg *= "The dimension should either be a Symbol or Integer indicating the direction "
    msg *= "of the condition!\n"
    msg *= "  x-direction: `:x` or `1`\n"
    msg *= "  y-direction: `:y` or `2`\n"
    msg *= "  z-direction: `:z` or `3`\n"
    return msg
end

function get_dim(dim::I) where {I<:Integer}
    in(dim, 1:3) || throw(ArgumentError(wrong_dim_err_msg(dim)))
    I <: UInt8 && return dim
    return convert(UInt8, dim)
end

function get_dim(dim::Symbol)
    if !haskey(SYMBOL_TO_DIM, dim)
        throw(ArgumentError(wrong_dim_err_msg(dim)))
    end
    return SYMBOL_TO_DIM[dim]
end

function get_dims(dimspec::Vector{T}) where {T<:Union{Integer,Symbol}}
    if length(dimspec) > 3
        throw(ArgumentError("too many dimensions specified!"))
    end
    dims::Vector{UInt8} = [get_dim(dim) for dim in dimspec]
    return dims
end

function add_boundary_condition!(body::AbstractBody, conditions::Vector{BC},
                                 condition::BC) where {BC<:AbstractCondition}
    check_boundary_condition_conflicts(body, condition)
    push!(conditions, condition)
    return nothing
end

function check_boundary_condition_conflicts(body::AbstractBody, condition::BC) where {BC}
    if has_boundary_condition_conflict(body, condition)
        msg = "the specified condition conflicts with already existing conditions!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function has_boundary_condition_conflict(body::AbstractBody, condition::BC) where {BC}
    for existing_condition in body.single_dim_bcs
        conditions_conflict(body, existing_condition, condition) && return true
    end
    for existing_condition in body.posdep_single_dim_bcs
        conditions_conflict(body, existing_condition, condition) && return true
    end
    for existing_condition in body.pos_single_dim_bcs
        conditions_conflict(body, existing_condition, condition) && return true
    end
    for existing_condition in body.data_bcs
        conditions_conflict(body, existing_condition, condition) && return true
    end
    return false
end

function conditions_conflict(body::AbstractBody, a::AbstractCondition, b::AbstractCondition)
    same_field = conditions_have_same_field(a, b)
    same_dim = conditions_have_same_dim(a, b)
    points_intersect = point_sets_intersect(body.point_sets, a.point_set, b.point_set)
    return same_field && same_dim && points_intersect
end

function conditions_have_same_field(a::AbstractCondition, b::AbstractCondition)
    return a.field === b.field
end

function conditions_have_same_dim(a::AbstractCondition, b::AbstractCondition)
    return a.dim == b.dim
end

function conditions_have_same_dim(a::DataBC, b::AbstractCondition)
    return b.dim in a.dims
end

function conditions_have_same_dim(a::AbstractCondition, b::DataBC)
    return a.dim in b.dims
end

function conditions_have_same_dim(a::DataBC, b::DataBC)
    return intersect(a.dims, b.dims) != UInt8[]
end

function check_boundary_condition_function(f::F) where {F<:Function}
    func_method = get_method_of_function(f)
    args = get_argument_names_of_function(func_method)
    if length(args) == 1 && args[1] === :t
        type = :sdbc
    elseif length(args) == 1 && args[1] === :p
        type = :psdbc
    elseif length(args) == 2 && args[1] === :p && args[2] === :t
        type = :pdsdbc
    else
        msg = "wrong arguments or argument names for condition function!\n"
        msg *= "Boundary conditions support only three type of functions:\n"
        msg *= "  `f(p, t)`\n"
        msg *= "  `f(t)`\n"
        msg *= "  `f(p)`\n"
        msg *= "with `t` beeing the current time and `p` a the position vector [x, y, z] "
        msg *= "of each point in the point set!\n"
        throw(ArgumentError(msg))
    end
    return type
end

function field_to_name(field::Symbol)
    if field === :velocity || field === :velocity_half
        name = "velocity"
    elseif field === :b_ext
        name = "force density"
    else
        name = string(field)
    end
    return name
end

"""
    velocity_bc!(fun, body, set_name, dim)

Specify velocity boundary condition for points of the set `set_name` in `body`.
The value of the boundary condition is calculated with the function `fun` at every time
step.

# Arguments

- `fun::Function`: Condition function for the calculation of a value, should return a
    `Float64`. If the condition function returns a `NaN`, then this value is ignored, which
    can be used to turn conditions off after a specified period of time. This function
    accepts one ore two positional arguments and is aware of the argument names.
    Possible arguments and names:
    - `fun(t)`: The function will receive the current time `t` at every time step.
        This makes it possible to specify conditions that change over time.
    - `fun(p, t)`: This function will be processed for every point of `set_name` and
        receives the reference position of a point as `SVector{3}` and the current time `t`
        at every time step. This makes it possible to specify conditions that
        also depend on the position of a point.
- `body::AbstractBody`: [`Body`](@ref) the condition is specified on.
- `set_name::Symbol`: The name of a point set of this body.
- `dim::Union{Integer,Symbol}`: Direction of the condition, either specified as Symbol or
    integer.
    -  x-direction: `:x` or `1`
    -  y-direction: `:y` or `2`
    -  z-direction: `:z` or `3`

# Throws

- Error if the body does not contain a set with `set_name`.
- Error if the direction is not correctly specified.
- Error if function is not suitable as condition function and has the wrong arguments.

# Example

```julia-repl
julia> velocity_bc!(t -> 2.0, body, :all_points, 1)

julia> velocity_bc!((p,t) -> p[1] * t, body, :all_points, :y)

julia> velocity_bc!(t -> t > 0.00001 ? 1.0 : NaN, body, :all_points, :z)

julia> body
1000-point Body{BBMaterial{NoCorrection}}:
  1 point set(s):
    1000-point set `all_points`
  3 boundary condition(s):
    BC on velocity: point_set=all_points, dim=1
    BC on velocity: point_set=all_points, dim=3
    Pos.-dep. BC on velocity: point_set=all_points, dim=2
```
"""
function velocity_bc!(f::F, body::AbstractBody, point_set::Symbol,
                      dimension::Union{Integer,Symbol}) where {F<:Function}
    check_if_set_is_defined(body.point_sets, point_set)
    type = check_boundary_condition_function(f)
    dim = get_dim(dimension)
    if type === :sdbc
        sdbc = SingleDimBC(f, :velocity_half, point_set, dim)
        add_boundary_condition!(body, body.single_dim_bcs, sdbc)
    elseif type === :pdsdbc
        pdsdbc = PosDepSingleDimBC(f, :velocity_half, point_set, dim)
        add_boundary_condition!(body, body.posdep_single_dim_bcs, pdsdbc)
    end
    return nothing
end

"""
    velocity_databc!(body, data, set_name, dims)

$(internal_api_warning())
$(experimental_api_warning())

Specifies velocity boundary conditions for points of the set `set_name` in `body`.
The value of the boundary condition is assigned by reading the corresponding positions
in the matrix `data`. Multiple dimensions can be handled at once.

!!! warning "Compatibility feature with other packages"
    This method / feature is used for compatibility with other packages developing with
    `Peridynamics.jl`. It is likely to change in the future, since the functionality of
    updating the values of the matrix during the simulation is not yet implemented.
    Consequently, at this stage, it is only available as a private API to facilitate future
    modifications and ensure easier implementation of changes.

# Arguments

- `body::AbstractBody`: [`Body`](@ref) the condition is specified on.
- `data::Matrix`: A matrix of size `length(dims) x n_points` that contains the values of
    the boundary condition for each point in the body. But only the conditions of points
    contained in the set `set_name` are applied during the simulation! It should be noted,
    that the value of the `data` matrix is constant and currently cannot be updated during
    the simulation. The data matrix is not checked for `NaN` values, since this is handled
    in the `apply_bc!` function. If it contains `NaN` values, then these values are ignored.
- `set_name::Symbol`: The name of a point set of this body. The condition applies only to
    the points in this set, even if the data matrix contains values for all points in the
    body.
- `dims::Vector{Union{Integer,Symbol}}`: Vector containing the directions of the condition
    that should be applied, either specified as Symbol or integer.
    -  x-direction: `:x` or `1`
    -  y-direction: `:y` or `2`
    -  z-direction: `:z` or `3`
    It should not contain more than 3 elements, and the elements should be unique. The order
    of the elements does not matter, however it must match the values in the data matrix.
    So if the first column of the data matrix contains the values for the x-direction,
    then the first element of `dims` should be `1` or `:x`, and so on.

# Throws

- Errors if the body does not contain a set with `set_name`.
- Errors if the directions are not correctly specified.
- Errors if the dimensions of the data matrix are incorrect.
"""
function velocity_databc!(body::AbstractBody, data::Matrix, name::Symbol,
                          dimspec::Vector{T}) where {T<:Union{Integer,Symbol}}
    check_if_set_is_defined(body.point_sets, name)
    check_databc_data_dimensions(body, data, dimspec)
    dims = get_dims(dimspec)
    databc = DataBC(data, :velocity_half, name, dims)
    add_boundary_condition!(body, body.data_bcs, databc)
    return nothing
end

function check_databc_data_dimensions(body::AbstractBody, data::Matrix, dimspec::Vector)
    if size(data, 2) != n_points(body)
        msg = "the data matrix has a different number of columns than the\n"
        msg *= "number of points in the body!\n"
        throw(ArgumentError(msg))
    end
    if size(data, 1) != length(dimspec)
        msg = "the data matrix has a different number of rows than the elements\n"
        msg *= "in the dimension matrix!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

"""
    forcedensity_bc!(fun, body, set, dim)

Specify force density boundary condition for points of the set `set_name` in `body`.
The value of the boundary condition is calculated with the function `fun` at every time
step.

# Arguments

- `fun::Function`: Condition function for the calculation of a value, should return a
    `Float64`. If the condition function returns a `NaN`, then this value is ignored, which
    can be used to turn conditions off after a specified period of time. This function
    accepts one ore two positional arguments and is aware of the argument names.
    Possible arguments and names:
    - `fun(t)`: The function will receive the current time `t` at every time step.
        This makes it possible to specify conditions that change over time.
    - `fun(p, t)`: This function will be processed for every point of `set_name` and
        receives the reference position of a point as `SVector{3}` and the current time `t`
        at every time step. This makes it possible to specify conditions that
        also depend on the position of a point.
- `body::AbstractBody`: [`Body`](@ref) the condition is specified on.
- `set_name::Symbol`: The name of a point set of this body.
- `dim::Union{Integer,Symbol}`: Direction of the condition, either specified as Symbol or
    integer.
    -  x-direction: `:x` or `1`
    -  y-direction: `:y` or `2`
    -  z-direction: `:z` or `3`

# Throws

- Error if the body does not contain a set with `set_name`.
- Error if the direction is not correctly specified.
- Error if function is not suitable as condition function and has the wrong arguments.

# Example

```julia-repl
julia> forcedensity_bc!(t -> 8000.0, body, :all_points, :x)

julia> forcedensity_bc!((p,t) -> p[1] * t, body, :all_points, :y)

julia> forcedensity_bc!(t -> t > 0.00001 ? 8000.0 : NaN, body, :all_points, :z)

julia> body
1000-point Body{BBMaterial{NoCorrection}}:
  1 point set(s):
    1000-point set `all_points`
  3 boundary condition(s):
    BC on force density: point_set=all_points, dim=1
    BC on force density: point_set=all_points, dim=3
    Pos.-dep. BC on force density: point_set=all_points, dim=2
```
"""
function forcedensity_bc!(f::F, body::AbstractBody, point_set::Symbol,
                          dimension::Union{Integer,Symbol}) where {F<:Function}
    check_if_set_is_defined(body.point_sets, point_set)
    type = check_boundary_condition_function(f)
    dim = get_dim(dimension)
    if type === :sdbc
        sdbc = SingleDimBC(f, :b_ext, point_set, dim)
        add_boundary_condition!(body, body.single_dim_bcs, sdbc)
    elseif type === :psdbc
        psdbc = PosSingleDimBC(f, :b_ext, point_set, dim)
        add_boundary_condition!(body, body.pos_single_dim_bcs, psdbc)
    elseif type === :pdsdbc
        pdsdbc = PosDepSingleDimBC(f, :b_ext, point_set, dim)
        add_boundary_condition!(body, body.posdep_single_dim_bcs, pdsdbc)
    end
    return nothing
end

"""
    forcedensity_databc!(body, data, set_name, dims)

$(internal_api_warning())
$(experimental_api_warning())

Specifies forcedensity boundary conditions for points of the set `set_name` in `body`.
The value of the boundary condition is assigned by reading the corresponding positions
in the matrix `data`. Multiple dimensions can be handled at once.

!!! warning "Compatibility feature with other packages"
    This method / feature is used for compatibility with other packages developing with
    `Peridynamics.jl`. It is likely to change in the future, since the functionality of
    updating the values of the matrix during the simulation is not yet implemented.
    Consequently, at this stage, it is only available as a private API to facilitate future
    modifications and ensure easier implementation of changes.

# Arguments

- `body::AbstractBody`: [`Body`](@ref) the condition is specified on.
- `data::Matrix`: A matrix of size `length(dims) x n_points` that contains the values of
    the boundary condition for each point in the body. But only the conditions of points
    contained in the set `set_name` are applied during the simulation! It should be noted,
    that the value of the `data` matrix is constant and currently cannot be updated during
    the simulation. The data matrix is not checked for `NaN` values, since this is handled
    in the `apply_bc!` function. If it contains `NaN` values, then these values are ignored.
- `set_name::Symbol`: The name of a point set of this body. The condition applies only to
    the points in this set, even if the data matrix contains values for all points in the
    body.
- `dims::Vector{Union{Integer,Symbol}}`: Vector containing the directions of the condition
    that should be applied, either specified as Symbol or integer.
    -  x-direction: `:x` or `1`
    -  y-direction: `:y` or `2`
    -  z-direction: `:z` or `3`
    It should not contain more than 3 elements, and the elements should be unique. The order
    of the elements does not matter, however it must match the values in the data matrix.
    So if the first column of the data matrix contains the values for the x-direction,
    then the first element of `dims` should be `1` or `:x`, and so on.

# Throws

- Errors if the body does not contain a set with `set_name`.
- Errors if the directions are not correctly specified.
- Errors if the dimensions of the data matrix are incorrect.
"""
function forcedensity_databc!(body::AbstractBody, data::Matrix, name::Symbol,
                              dimspec::Vector{T}) where {T<:Union{Integer,Symbol}}
    check_if_set_is_defined(body.point_sets, name)
    check_databc_data_dimensions(body, data, dimspec)
    dims = get_dims(dimspec)
    databc = DataBC(data, :b_ext, name, dims)
    add_boundary_condition!(body, body.data_bcs, databc)
    return nothing
end



"""
    displacement_bc!(f::Function, body::Body, points::Vector{Int}, dim::Int)

$(internal_api_warning())
$(experimental_api_warning())

!!! warning "Compatibility limited"
    This boundary condition type only works with the NewtonKrylov solver for now.

Apply a displacement boundary condition to specified points in a given dimension.
A factor `β` is calculated as `n / n_steps`, where `n` is the current step number and
`n_steps` is the total number of steps. Then the return value of the function `f` is
multiplicated with `β` at each time step. This means that the displacement will be applied
gradually over the course of the simulation.

# Arguments:
- `fun::Function`: Condition function for the calculation of a value, should return a
    `Float64`. If the condition function returns a `NaN`, then this value is ignored, which
    can be used to turn conditions off. This function accepts only one argument and is aware
    of the argument names. Possible arguments and names:
    - `fun(p)`: This function will be processed for every point of `set_name` and
        receives the reference position of a point as `SVector{3}` at every time step.
        This makes it possible to specify conditions that depend on the position of a point.
- `body::AbstractBody`: [`Body`](@ref) the condition is specified on.
- `set_name::Symbol`: The name of a point set of this body.
- `dim::Union{Integer,Symbol}`: Direction of the condition, either specified as Symbol or
    integer.
    -  x-direction: `:x` or `1`
    -  y-direction: `:y` or `2`
    -  z-direction: `:z` or `3`

# Returns:
- `nothing`: No return value. The boundary condition is added to the body.

# Example:
```julia
# Apply constant prescribed displacement
points = [1, 2, 3, 4, 5]
displacement_bc!(p -> 0.1, body, points, 1)

# Apply prescribed displacement based on position
displacement_bc!(p -> 0.01 * p[1], body, points, 2)
```
"""
function displacement_bc!(f::F, body::AbstractBody, point_set::Symbol,
                          dimension::Union{Integer,Symbol}) where {F}
    check_if_set_is_defined(body.point_sets, point_set)
    type = check_boundary_condition_function(f)
    if type !== :psdbc
        msg = "the function for a displacement boundary condition must be of the form\n"
        msg *= "`f(p)`!\n"
        throw(ArgumentError(msg))
    end
    dim = get_dim(dimension)
    psdbc = PosSingleDimBC(f, :displacement, point_set, dim)
    add_boundary_condition!(body, body.pos_single_dim_bcs, psdbc)
    return nothing
end
