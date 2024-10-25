function get_correction(mat::AbstractBondSystemMaterial{<:AbstractZEMStabilization},
                        ::Int, ::Int, ::Int)
    return mat.zem_stabilization
end

"""
    ZEMSilling(; Cs=0.0)

TODO
"""
struct ZEMSilling <: AbstractZEMStabilization
    Cs::Float64
    function ZEMSilling(; Cs::Real=100.0)
        return new(Cs)
    end
end
