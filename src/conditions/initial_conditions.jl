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

function apply_ic!(b::AbstractBodyChunk, ic::SingleDimIC)
    for point_id in b.psets[ic.point_set]
        setindex!(getfield(b.store, ic.field), ic.value, ic.dim, point_id)
    end
    return nothing
end
