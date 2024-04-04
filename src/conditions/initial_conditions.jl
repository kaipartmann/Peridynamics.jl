"""

"""
struct SingleDimIC <: AbstractCondition
    value::Float64
    field::Symbol
    point_set::Symbol
    dim::UInt8
end

function override_eachother(a::SingleDimIC, b::SingleDimIC)
    same_field = a.field === b.field
    same_point_set = a.point_set === b.point_set
    same_dim = a.dim == b.dim
    return same_field && same_point_set && same_dim
end

function apply_initial_conditions!(b::AbstractBodyChunk, body::AbstractBody)
    apply_single_dim_ic!(b, body)
    return nothing
end

@inline function apply_single_dim_ic!(b::AbstractBodyChunk, body::AbstractBody)
    for ic in body.single_dim_ics
        apply_ic!(b, ic)
    end
    return nothing
end

function apply_ic!(b::AbstractBodyChunk, ic::SingleDimIC)
    for point_id in b.psets[ic.point_set]
        setindex!(get_storage_field(b.storage, ic.field), ic.value, ic.dim, point_id)
    end
    return nothing
end

"""
    velocity_ic!(b::AbstractBody, name::Symbol, d::DimensionSpec, value::Real)

specifies initital conditions for the velocity of points in point set `name` on body `b`

# Arguments

- `b::AbstractBody`: peridynamic body
- `name::Symbol`: point set on body `b`
- `d::DimensionSpec`: direction of velocity
- `value::Real`: initial velocity value

# Throws

- error if no point set called `name` exists
- error if dimension is not correctly specified

# Example

```julia-repl
julia> velocity_ic!(b, :set_b, :y, 20)
julia> b.single_dim_ics
1-element Vector{Peridynamics.SingleDimIC}:
 Peridynamics.SingleDimIC(20.0, :velocity, :set_b, 0x02)
```
"""
function velocity_ic!(b::AbstractBody, name::Symbol, d::Union{Integer,Symbol}, value::Real)
    check_if_set_is_defined(b.point_sets, name)
    dim = get_dim(d)
    sdic = SingleDimIC(convert(Float64, value), :velocity, name, dim)
    _condition!(b.single_dim_ics, sdic)
    return nothing
end
