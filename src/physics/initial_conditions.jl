"""

"""
struct SingleDimIC <: AbstractCondition
    value::Float64
    field::Symbol
    point_set::Symbol
    dim::UInt8
end

function Base.show(io::IO, bc::SingleDimIC)
    print(io, "SingleDimIC")
    print(io, msg_fields_in_brackets(bc, (:value, :field, :point_set, :dim)))
    return nothing
end

function override_eachother(a::SingleDimIC, b::SingleDimIC)
    same_field = a.field === b.field
    same_point_set = a.point_set === b.point_set
    same_dim = a.dim == b.dim
    return same_field && same_point_set && same_dim
end

function apply_initial_conditions!(b::AbstractBodyChunk, body::AbstractBody)
    for ic in body.single_dim_ics
        apply_ic!(b, ic)
    end
    for ic in body.posdep_single_dim_ics
        apply_ic!(b, ic)
    end
    return nothing
end

function apply_ic!(b::AbstractBodyChunk, ic::SingleDimIC)
    for point_id in b.psets[ic.point_set]
        setindex!(get_point_data(b.storage, ic.field), ic.value, ic.dim, point_id)
    end
    return nothing
end

struct PosDepSingleDimIC{F<:Function} <: AbstractCondition
    fun::F
    field::Symbol
    point_set::Symbol
    dim::UInt8
end

function Base.show(io::IO, bc::PosDepSingleDimIC)
    print(io, "PosDepSingleDimIC")
    print(io, msg_fields_in_brackets(bc, (:field, :point_set, :dim)))
    return nothing
end

@inline function (b::PosDepSingleDimIC{F})(p::AbstractVector) where {F}
    value::Float64 = b.fun(p)
    return value
end

function override_eachother(a::PosDepSingleDimIC, b::PosDepSingleDimIC)
    same_field = a.field === b.field
    same_point_set = a.point_set === b.point_set
    same_dim = a.dim == b.dim
    return same_field && same_point_set && same_dim
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

"""
    velocity_ic!(body, set, dim, value)
    velocity_ic!(fun, body, set, dim)

Specifies initital conditions for the velocity of points in point set `set` on `body`

# Arguments

- `body::AbstractBody`: Peridynamic body
- `set::Symbol`: Point set on `body`
- `dim::Union{Integer,Symbol}`: Direction of velocity
- `value::Real`: Initial velocity value
- `fun::Function`: Velocity condition function #TODO

# Throws

- Error if no point set called `set` exists
- Error if dimension is not correctly specified

# Example

```julia-repl
julia> velocity_ic!(b, :set_b, :y, 20)

julia> b.single_dim_ics
1-element Vector{Peridynamics.SingleDimIC}:
 Peridynamics.SingleDimIC(20.0, :velocity, :set_b, 0x02)
```
"""
function velocity_ic! end

function velocity_ic!(body::AbstractBody, point_set::Symbol,
                      dimension::Union{Integer,Symbol}, value::Real)
    check_if_set_is_defined(body.point_sets, point_set)
    sdic = SingleDimIC(convert(Float64, value), :velocity, point_set, get_dim(dimension))
    add_condition!(body.single_dim_ics, sdic)
    return nothing
end

function velocity_ic!(f::F, body::AbstractBody, point_set::Symbol,
                      dimension::Union{Integer,Symbol}) where {F<:Function}
    check_if_set_is_defined(body.point_sets, point_set)
    check_initial_condition_function(f)
    sdic = PosDepSingleDimIC(f, :velocity, point_set, get_dim(dimension))
    add_condition!(body.posdep_single_dim_ics, sdic)
    return nothing
end

function check_initial_condition_function(f::F) where {F<:Function}
    func_method = get_method_of_function(f)
    args = get_argument_names_of_function(func_method)
    if length(args) != 1 || args[1] !== :p
        msg = "wrong arguments for position dependent initial condition function!\n"
        msg *= "Initial condition functions support only `p` as argument name for the\n"
        msg *= "position of a point in the point set.\n"
        throw(ArgumentError(msg))
    end
    return nothing
end
