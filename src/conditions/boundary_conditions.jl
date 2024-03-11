
"""

"""
struct SingleDimBC{F<:Function,G<:Function} <: AbstractCondition
    fun::F
    cff::G
    field::Symbol
    point_set::Symbol
    dim::UInt8
end

@inline function (b::SingleDimBC{F,G})(t::Float64) where {F<:Function,G<:Function}
    value::Float64 = b.fun(t)
    return value
end

@inline function (b::SingleDimBC{F,G})(s::AbstractStorage) where {F<:Function,G<:Function}
    return b.cff(s)
end

@inline get_cff_velocity_half(s::AbstractStorage) = s.velocity_half
@inline get_cff_velocity(s::AbstractStorage) = s.velocity
@inline get_cff_b_ext(s::AbstractStorage) = s.b_ext

function override_eachother(a::SingleDimBC, b::SingleDimBC)
    same_field = a.field === b.field
    same_point_set = a.point_set === b.point_set
    same_dim = a.dim == b.dim
    return same_field && same_point_set && same_dim
end

function apply_bcs!(b::AbstractBodyChunk, time::Float64)
    for bc in b.sdbcs
        apply_bc!(b.store, b.psets, bc, time)
    end
    return nothing
end

function apply_bc!(s::AbstractStorage, psets::Dict{Symbol,Vector{Int}},
                   bc::SingleDimBC{F,G}, time::Float64) where {F,G}
    value = bc(time)
    isnan(value) && return nothing
    apply_sdbc!(bc(s), value, bc.dim, psets[bc.point_set])
    return nothing
end

@inline function apply_sdbc!(field::Matrix{Float64}, value::Float64, dim::UInt8,
                             point_ids::Vector{Int})
    @simd for i in point_ids
        @inbounds field[dim, i] = value
    end
    return nothing
end
