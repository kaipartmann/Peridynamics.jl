
"""

"""
struct SingleDimBC{F<:Function} <: AbstractCondition
    fun::F
    field::Symbol
    point_set::Symbol
    dim::UInt8
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

function apply_bc!(s::AbstractStorage, psets::Dict{Symbol,Vector{Int}}, bc::SingleDimBC,
                   time::Float64)
    value = bc.fun(time)
    isnan(value) && return nothing
    for point_id in psets[bc.point_set]
        setindex!(getfield(s, bc.field), value, bc.dim, point_id)
    end
    return nothing
end
