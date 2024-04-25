struct SingleDimBC{F<:Function} <: AbstractCondition
    fun::F
    field::Symbol
    point_set::Symbol
    dim::UInt8
end

@inline function (b::SingleDimBC{F})(t::Float64) where {F}
    value::Float64 = b.fun(t)
    return value
end

function override_eachother(a::SingleDimBC, b::SingleDimBC)
    same_field = a.field === b.field
    same_point_set = a.point_set === b.point_set
    same_dim = a.dim == b.dim
    return same_field && same_point_set && same_dim
end

function apply_bcs!(b::AbstractBodyChunk, time::Float64)
    for bc in b.sdbcs
        apply_bc!(b.storage, b.psets, bc, time)
    end
    for bc in b.pdsdbcs
        apply_bc!(b.storage, b.psets, bc, b.system.position, time)
    end
    return nothing
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

@inline function (b::PosDepSingleDimBC{F})(p::AbstractVector, t::Float64) where {F}
    value::Float64 = b.fun(p, t)
    return value
end

function override_eachother(a::PosDepSingleDimBC, b::PosDepSingleDimBC)
    same_field = a.field === b.field
    same_point_set = a.point_set === b.point_set
    same_dim = a.dim == b.dim
    return same_field && same_point_set && same_dim
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
        value = bc(SVector{3}(pos[1,i], pos[2,i], pos[3,i]), t)
        if !isnan(value)
            field[bc.dim, i] = value
        end
    end
    return nothing
end

function get_dim(dim::I) where {I<:Integer}
    in(dim, 1:3) || error("specified dimension should be 1=x, 2=y, or 3=z!\n")
    isa(UInt8, I) && return dim
    return convert(UInt8, dim)
end

function get_dim(dim::Symbol)
    if !haskey(SYMBOL_TO_DIM, dim)
        error("unknown dimension symbol $(dim)! Should be :x, :y, or :z!\n")
    end
    return SYMBOL_TO_DIM[dim]
end

function _condition!(conditions::Vector{B}, condition::B) where {B<:AbstractCondition}
    # check if conditions override each other!
    if is_duplicate(condition, conditions)
        error("duplicate conditions for point set $(condition.point_set)!\n")
    end
    push!(conditions, condition)
    return nothing
end

"""
    velocity_bc!(fun, body, set, dim)

Specifies velocity boundary conditions for point set `set` on `body`

# Arguments

- `fun::Function`: Velocity condition function
- `body::AbstractBody`: Peridynamic body
- `set::Symbol`: Point set on `body`
- `dim::DimensionSpec`: Direction of velocity

# Throws

- Error if no point set called `set` exists
- Error if dimension is not correctly specified
- Error if function is not suitable as condition function

# Example

```julia-repl
julia> velocity_bc!(t -> -9.81 * t, b, :set_bottom, :y)

julia> b.single_dim_bcs

1-element Vector{Peridynamics.SingleDimBC}:
 Peridynamics.SingleDimBC{var"#15#16"}(var"#15#16"(), :velocity_half, :set_bottom, 0x02)

julia> velocity_bc!(t -> 40, b, :set_a, 1)

julia> b.single_dim_bcs
2-element Vector{Peridynamics.SingleDimBC}:
 Peridynamics.SingleDimBC{var"#15#16"}(var"#15#16"(), :velocity_half, :set_bottom, 0x02)
 Peridynamics.SingleDimBC{var"#17#18"}(var"#17#18"(), :velocity_half, :set_a, 0x01)
```
"""
function velocity_bc!(f::F, b::AbstractBody, name::Symbol, d::Union{Integer,Symbol}) where {F<:Function}
    check_if_set_is_defined(b.point_sets, name)
    type = check_condition_function(f)
    dim = get_dim(d)
    if type === :sdbc
        sdbc = SingleDimBC(f, :velocity_half, name, dim)
        _condition!(b.single_dim_bcs, sdbc)
    elseif type === :pdsdbc
        pdsdbc = PosDepSingleDimBC(f, :velocity_half, name, dim)
        _condition!(b.posdep_single_dim_bcs, pdsdbc)
    end
    return nothing
end

"""
    forcedensity_bc!(fun, body, set, dim)

Specifies boundary conditions for force density on points of point set `set` on `body`

# Arguments

- `fun::Function`: Condition function
- `body::AbstractBody`: Peridynamic body
- `set::Symbol`: Point set on `body`
- `dim::DimensionSpec`: Direction of force density

# Throws

- Error if no point set called `set` exists
- Error if dimension is not correctly specified
- Error if function is not suitable as condition function

# Example

```julia-repl
julia> forcedensity_bc!(t -> 40, b, :set_a, 1)

julia> b.single_dim_bcs
1-element Vector{Peridynamics.SingleDimBC}:
 Peridynamics.SingleDimBC{var"#25#26"}(var"#25#26"(), :b_ext, :set_a, 0x01)
```
"""
function forcedensity_bc!(f::F, b::AbstractBody, name::Symbol, d::Union{Integer,Symbol}) where {F<:Function}
    check_if_set_is_defined(b.point_sets, name)
    type = check_condition_function(f)
    dim = get_dim(d)
    if type === :sdbc
        sdbc = SingleDimBC(f, :b_ext, name, dim)
        _condition!(b.single_dim_bcs, sdbc)
    elseif type === :pdsdbc
        pdsdbc = PosDepSingleDimBC(f, :b_ext, name, dim)
        _condition!(b.posdep_single_dim_bcs, pdsdbc)
    end
    return nothing
end
