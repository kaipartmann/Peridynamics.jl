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
    isnan(value) && return nothing
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
        if !isnan(value)
            field[bc.dim, i] = value
        end
    end
    return nothing
end

struct VBC <: AbstractCondition
    mbc::AbstractMatrix
    field::Symbol
    point_set::Symbol
    dims::Vector{UInt8}
end

function Base.show(io::IO, @nospecialize(bc::VBC))
    print(io, "BC on ", field_to_name(bc.field), ": ")
    print(io, msg_fields_inline(bc, (:point_set)))
    return nothing
end

function apply_bc!(s::AbstractStorage, psets::Dict{Symbol,Vector{Int}}, bc::VBC, all_point_ids::Vector{Int})
    apply_sdbc!(get_point_data(s, bc.field), bc, psets[bc.point_set], all_point_ids)
    return nothing
end

@inline function apply_sdbc!(field::Matrix{Float64}, bc::VBC, point_ids::Vector{Int}, all_point_ids::Vector{Int})
    @simd for i in point_ids
        for j in eachindex(bc.dims)
            value = bc.mbc[j, all_point_ids[i]]
            dim = bc.dims[j]
            field[dim, i] = value
        end
    end
    return nothing
end

function apply_boundary_conditions!(b::AbstractBodyChunk, time::Float64)
    for bc in b.sdbcs
        apply_bc!(b.storage, b.psets, bc, time)
    end
    for bc in b.pdsdbcs
        apply_bc!(b.storage, b.psets, bc, b.system.position, time)
    end
    for bc in b.vbcs
        apply_bc!(b.storage, b.psets, bc, b.system.chunk_handler.point_ids)
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
    return false
end

function conditions_conflict(body::AbstractBody, a::AbstractCondition, b::AbstractCondition)
    same_field = a.field === b.field
    same_dim = a.dim == b.dim
    points_intersect = point_sets_intersect(body.point_sets, a.point_set, b.point_set)
    return same_field && same_dim && points_intersect
end

function check_boundary_condition_function(f::F) where {F<:Function}
    func_method = get_method_of_function(f)
    args = get_argument_names_of_function(func_method)
    if length(args) == 1 && args[1] === :t
        type = :sdbc
    elseif length(args) == 2 && args[1] === :p && args[2] === :t
        type = :pdsdbc
    else
        msg = "wrong arguments or argument names for condition function!\n"
        msg *= "Boundary conditions support only two type of functions:\n"
        msg *= "  `f(p, t)`\n"
        msg *= "  `f(t)`\n"
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

Specifies velocity boundary conditions for points of the set `set_name` in `body`.
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

- Errors if the body does not contain a set with `set_name`.
- Errors if the direction is not correctly specified.
- Errors if function is not suitable as condition function and has the wrong arguments.

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

# A new method to apply velocity boundary conditions

velocity_bc!(v_matrix, body, set_name, dims)

Specifies velocity boundary conditions for points of the set `set_name` in `body`.

The value of the boundary condition is assigned by reading the corresponding positions
    in the matrix(`v_matrix`). 

Multiple dimensions can be handled at once.

More importantly, if we cannot describe the evolution of boundary conditions over time
    using a function, and can only obtain the boundary values directly 
    (e.g., in fluid-structure interaction problems), 
    or if the boundary changes—such as when the original boundary disappears and a new one forms 
    (e.g., in cases of surface ablation or melting)—we can simply update the corresponding values
    in the `v_matrix` at each iteration to implement this.

Of course, this method of application is more suitable for pressure boundaries
    from a physical perspective.  

# Example
julia> dims = UInt8[1, 2, 3] 
julia> v_matrix = zeros(3, body.n_points)
julia> velocity_bc!(v_matrix, body, :all_points, dims)

Here, the elements of dims are of type UInt8.     
    -  x-direction: `1`
    -  y-direction: `2`
    -  z-direction: `3`

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

function velocity_bc!(f::Matrix, b::AbstractBody, name::Symbol, dims::Vector{UInt8})
    vbc = VBC(f, :velocity_half, name, dims)
    add_boundary_condition!(b, b.v_bcs, vbc)
    return nothing
end

"""
    forcedensity_bc!(fun, body, set, dim)

Specifies force density boundary conditions for points of the set `set_name` in `body`.
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

- Errors if the body does not contain a set with `set_name`.
- Errors if the direction is not correctly specified.
- Errors if function is not suitable as condition function and has the wrong arguments.

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

# A new method to apply forcedensity boundary conditions

forcedensity_bc!(f_matrix, body, set_name, dims)

Specifies forcedensity boundary conditions for points of the set `set_name` in `body`.

The value of the boundary condition is assigned by reading the corresponding positions
    in the matrix(`f_matrix`). 

Multiple dimensions can be handled at once.

More importantly, if we cannot describe the evolution of boundary conditions over time
    using a function, and can only obtain the boundary values directly 
    (e.g., in fluid-structure interaction problems), 
    or if the boundary changes—such as when the original boundary disappears and a new one forms 
    (e.g., in cases of surface ablation or melting)—we can simply update the corresponding values
    in the `v_matrix` at each iteration to implement this. 

By applying the boundary conditions through the matrix, we can adjust the conditions for each PD node
    at every step, offering greater flexibility and control.

# Example
julia> dims = UInt8[1, 2, 3] 
julia> f_matrix = zeros(3, body.n_points)
julia> forcedensity_bc!(f_matrix, body, :all_points, dims)

Here, the elements of dims are of type UInt8.     
    -  x-direction: `1`
    -  y-direction: `2`
    -  z-direction: `3`

"""
function forcedensity_bc!(f::F, body::AbstractBody, point_set::Symbol,
                          dimension::Union{Integer,Symbol}) where {F<:Function}
    check_if_set_is_defined(body.point_sets, point_set)
    type = check_boundary_condition_function(f)
    dim = get_dim(dimension)
    if type === :sdbc
        sdbc = SingleDimBC(f, :b_ext, point_set, dim)
        add_boundary_condition!(body, body.single_dim_bcs, sdbc)
    elseif type === :pdsdbc
        pdsdbc = PosDepSingleDimBC(f, :b_ext, point_set, dim)
        add_boundary_condition!(body, body.posdep_single_dim_bcs, pdsdbc)
    end
    return nothing
end

function forcedensity_bc!(f::Matrix, b::AbstractBody, name::Symbol, dims::Vector{UInt8})
    vbc = VBC(f, :b_ext, name, dims)
    add_boundary_condition!(b, b.v_bcs, vbc)
    return nothing
end
