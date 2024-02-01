

"""
    SVBCSetup(field, point_set, dim, values_in_time)

The setup for a single value condition.

# Fields
- `field::Symbol`: field in storage that will recieve the value
- `point_set::Vector{Int}`: point id's of the field in storage that will recieve the value
- `dim::Int`: dimension (1=x, 2=y, 3=z) of the field in storage that will recieve the value
"""
struct SingleValueBCSetup
    field::Symbol
    point_set::Vector{Int}
    dim::Int
    values_in_time::Vector{Float64}
end

function apply_bc!(s::S, svbc::SingleValueBCSetup, n::Int) where {S<:AbstractStorage}
    value = svbc.values_in_time[n] # n = time step
    isnan(value) && return nothing
    field = getfield(s, svbc.field)
    for i in svbc.point_set
        setindex!(field, value, svbc.dim, i)
    end
    return nothing
end

function check_single_value_bc_func(f::F) where {F<:Function}
    func_method = get_method_of_function(f)

end

struct VelocityBC{F<:Function} <: AbstractSingleValueBC
    fun::F
    point_set::Symbol
    dim::UInt8
end

function get_bc_setup(bc::VelocityBC, times::V) where {V<:AbstractVector{<:Real}}
    values_in_time = bc.fun.(times)
    return SingleValueBCSetup(:velocity_half, bc.point_set, bc.dim, values_in_time)
end

struct ForceDensityBC{F<:Function} <: AbstractSingleValueBC
    fun::F
    point_set::Symbol
    dim::UInt8
end

function get_bc_setup(bc::ForceDensityBC, times::V) where {V<:AbstractVector{<:Real}}
    values_in_time = bc.fun.(times)
    return SingleValueBCSetup(:b_ext, bc.point_set, bc.dim, values_in_time)
end
