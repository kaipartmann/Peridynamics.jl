function get_correction(mat::AbstractBondSystemMaterial{<:AbstractZEMStabilization},
                        ::Int, ::Int, ::Int)
    return mat.zem_stabilization
end

"""
    ZEMSilling(; Cs)

Zero-energy mode stabilization algorithm of Silling (2017). This is necessary for the
correspondence formulation to stabilize the zero-energy modes.
See also [`CMaterial`](@ref) on how to use this stabilization algorithm.

# Keywords
- `Cs::Real`: Stabilization factor. (default: `100.0`)
"""
struct ZEMSilling <: AbstractZEMStabilization
    Cs::Float64
    function ZEMSilling(; Cs::Real=100.0)
        return new(Cs)
    end
end
