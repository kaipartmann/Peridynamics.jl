struct SingleDimIC <: AbstractCondition
    value::Float64
    field::Symbol
    point_set::Symbol
    dim::UInt8
end

function Base.show(io::IO, @nospecialize(ic::SingleDimIC))
    print(io, "IC on ", field_to_name(ic.field), ": ")
    print(io, msg_fields_inline(ic, (:point_set, :dim)))
    return nothing
end

function apply_ic!(chunk::AbstractBodyChunk, ic::SingleDimIC)
    for point_id in chunk.psets[ic.point_set]
        setindex!(get_point_data(chunk.storage, ic.field), ic.value, ic.dim, point_id)
    end
    return nothing
end

struct PosDepSingleDimIC{F<:Function} <: AbstractCondition
    fun::F
    field::Symbol
    point_set::Symbol
    dim::UInt8
end

function Base.show(io::IO, @nospecialize(ic::PosDepSingleDimIC))
    print(io, "Pos.-dep. IC on ", field_to_name(ic.field), ": ")
    print(io, msg_fields_inline(ic, (:point_set, :dim)))
    return nothing
end

@inline function (ic::PosDepSingleDimIC{F})(p::AbstractVector) where {F}
    value::Float64 = ic.fun(p)
    return value
end

function apply_ic!(chunk::AbstractBodyChunk, ic::PosDepSingleDimIC)
    (; system, storage, psets) = chunk
    field = get_point_data(storage, ic.field)
    @simd for point_id in psets[ic.point_set]
        value = ic(get_coordinates(system, point_id))
        if !isnan(value)
            @inbounds setindex!(field, value, ic.dim, point_id)
        end
    end
    return nothing
end

function apply_initial_conditions!(chunk::AbstractBodyChunk, body::AbstractBody)
    for ic in body.single_dim_ics
        apply_ic!(chunk, ic)
    end
    for ic in body.posdep_single_dim_ics
        apply_ic!(chunk, ic)
    end
    return nothing
end

function check_initial_condition_function(f::F) where {F<:Function}
    func_method = get_method_of_function(f)
    args = get_argument_names_of_function(func_method)
    if length(args) != 1 || args[1] !== :p
        msg = "wrong arguments for position dependent initial condition function!\n"
        msg *= "Initial conditions support only functions of type `f(p)`, where "
        msg *= "`p` is the position vector [x, y, z] of each point in the point set!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function add_initial_condition!(body::AbstractBody, conditions::Vector{BC},
                                condition::BC) where {BC<:AbstractCondition}
    check_initial_condition_conflicts(body, condition)
    push!(conditions, condition)
    return nothing
end

function check_initial_condition_conflicts(body::AbstractBody, condition::BC) where {BC}
    if has_initial_condition_conflict(body, condition)
        msg = "the specified condition conflicts with already existing conditions!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function has_initial_condition_conflict(body::AbstractBody, condition::BC) where {BC}
    for existing_condition in body.single_dim_ics
        conditions_conflict(body, existing_condition, condition) && return true
    end
    for existing_condition in body.posdep_single_dim_ics
        conditions_conflict(body, existing_condition, condition) && return true
    end
    return false
end

"""
    velocity_ic!(body, set_name, dim, value)
    velocity_ic!(fun, body, set_name, dim)

Specifies velocity initial conditions for points of the set `set_name` in `body`.
The `value` of the initial condition is specified before time integration.
If a function `fun` is specified, then the value is with that function.

# Arguments

- `body::AbstractBody`: [`Body`](@ref) the condition is specified on.
- `set_name::Symbol`: The name of a point set of this body.
- `dim::Union{Integer,Symbol}`: Direction of the condition, either specified as Symbol or
    integer.
    -  x-direction: `:x` or `1`
    -  y-direction: `:y` or `2`
    -  z-direction: `:z` or `3`
- `value::Real`: Value that is specified before time integration.
- `fun::Function`: Condition function for the calculation of a value, should return a
    `Float64`. If the condition function returns a `NaN`, then this value is ignored, which
    can be used to turn off the condition for a specified position. This function
    accepts one ore two positional arguments and is aware of the argument names.
    Possible arguments and names:
    - `fun(p)`: The function will receive the reference position `p` of a point as
        `SVector{3}`.

# Throws

- Errors if the body does not contain a set with `set_name`.
- Errors if the direction is not correctly specified.
- Errors if function is not suitable as condition function and has the wrong arguments.

# Example

```julia-repl
julia> velocity_ic!(body, :all_points, :x, -100.0)

julia> body
1000-point Body{BBMaterial{NoCorrection}}:
  1 point set(s):
    1000-point set `all_points`
  1 initial condition(s):
    IC on velocity: point_set=all_points, dim=1
```
"""
function velocity_ic! end

function velocity_ic!(body::AbstractBody, point_set::Symbol,
                      dimension::Union{Integer,Symbol}, value::Real)
    check_if_set_is_defined(body.point_sets, point_set)
    sdic = SingleDimIC(convert(Float64, value), :velocity, point_set, get_dim(dimension))
    add_initial_condition!(body, body.single_dim_ics, sdic)
    return nothing
end

function velocity_ic!(f::F, body::AbstractBody, point_set::Symbol,
                      dimension::Union{Integer,Symbol}) where {F<:Function}
    check_if_set_is_defined(body.point_sets, point_set)
    check_initial_condition_function(f)
    pdsdic = PosDepSingleDimIC(f, :velocity, point_set, get_dim(dimension))
    add_initial_condition!(body, body.posdep_single_dim_ics, pdsdic)
    return nothing
end
