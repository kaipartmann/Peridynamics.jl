
"""

"""
struct SingleDimBC{F<:Function} <: AbstractCondition
    fun::F
    field::Symbol
    point_set::Symbol
    dim::UInt8
end

function (b::SingleDimBC{F})(t::Float64) where {F<:Function}
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
        apply_bc!(b.store, b.psets, bc, time)
    end
    return nothing
end

function apply_bc!(s::AbstractStorage, psets::Dict{Symbol,Vector{Int}}, bc::SingleDimBC{F},
                   time::Float64) where {F<:Function}
    value = bc(time)
    isnan(value) && return nothing
    apply_sdbc!(get_bc_field(s, bc.field), value, bc.dim, psets[bc.point_set])
    return nothing
end

@inline function apply_sdbc!(field::Matrix{Float64}, value::Float64, dim::UInt8,
                             point_ids::Vector{Int})
    @simd for i in point_ids
        @inbounds field[dim, i] = value
    end
    return nothing
end

@inline function get_bc_field(s::AbstractStorage, fieldname::Symbol)
    return getfield(s, fieldname)::Matrix{Float64}
end
