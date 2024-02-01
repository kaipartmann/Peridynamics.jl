
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
