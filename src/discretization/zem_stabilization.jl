function get_correction(mat::AbstractBondSystemMaterial{<:AbstractZEMStabilization},
                        ::Int, ::Int, ::Int)
    return mat.zem_stabilization
end

"""
    ZEMSilling(; Cs=0.0)

TODO
"""
@kwdef struct ZEMSilling <: AbstractZEMStabilization
    Cs::Float64 = 0.0
end
