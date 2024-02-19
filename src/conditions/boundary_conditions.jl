
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

function apply_bc!(b::AbstractBodyChunk, bc::SingleDimBC, time::Float64)
    value = bc.fun(time)
    isnan(value) && return nothing
    for point_id in b.point_sets[bc.point_set]
        setindex!(getfield(b.storage, bc.field), value, bc.dim, point_id)
    end
    return nothing
end
