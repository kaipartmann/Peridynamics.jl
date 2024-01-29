
struct VelocityIC{V<:AbstractVector{<:Int}} <: AbstractSingleValueIC
    value::Float64
    point_set::V
    dim::Int
end

function apply_ic!(s::S, ic::VelocityIC) where {S<:AbstractStorage}
    value = ic.value
    dim = ic.dim
    for i in ic.point_set
        s.velocity[dim, i] = value
    end
    return nothing
end

struct ForceDensityIC <: AbstractSingleValueIC
    value::Float64
    point_set::Vector{Int}
    dim::Int
end

function apply_ic!(s::S, ic::ForceDensityIC) where {S<:AbstractStorage}
    value = ic.value
    dim = ic.dim
    for i in ic.point_set
        s.b_ext[dim, i] = value
    end
    return nothing
end
