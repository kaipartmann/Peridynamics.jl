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

function apply_initial_conditions!(b::AbstractBodyChunk, body::Body)
    apply_single_dim_ic!(b, body)
    return nothing
end

@inline function apply_single_dim_ic!(b::AbstractBodyChunk, body::Body)
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
