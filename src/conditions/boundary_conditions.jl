
"""

"""
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
        apply_bc!(b.store, b.psets, bc, time)
    end
    for bc in b.pdsdbcs
        apply_bc!(b.store, b.psets, bc, b.system.position, time)
    end
    return nothing
end

function apply_bc!(s::AbstractStorage, psets::Dict{Symbol,Vector{Int}},
                   bc::SingleDimBC{F}, time::Float64) where {F}
    value = bc(time)
    isnan(value) && return nothing
    apply_sdbc!(get_storage_field(s, bc.field), value, bc.dim, psets[bc.point_set])
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
    apply_pdsdbc!(get_storage_field(s, bc.field), position, bc, psets[bc.point_set], time)
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
