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

"""
    ZEMWan()

Zero-energy mode stabilization algorithm of Wan et al. (2019), which is an improvement to
Silling's algorithm that does not require a stabilization parameter.
See also [`CMaterial`](@ref) on how to use this stabilization algorithm.
"""
struct ZEMWan <: AbstractZEMStabilization end

function calc_zem_stiffness_tensor(C, Kinv)
    C_1 = zero(SMatrix{3,3,Float64,9})
    for i in 1:3, j in 1:3
        sum_value = 0.0
        for k in 1:3, l in 1:3
            @inbounds sum_value += C[i, j, k, l] * Kinv[k, l]
        end
        C_1 = setindex(C_1, sum_value, i, j)
    end
    return C_1
end

function calc_rotated_zem_stiffness_tensor!(C_rotated, C, Kinv, R)
    for m in 1:3, n in 1:3, o in 1:3, p in 1:3
        sum_value = 0.0
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            @inbounds sum_value += C[i, j, k, l] * R[i, m] * R[j, n] * R[k, o] * R[l, p]
        end
        @inbounds C_rotated[m, n, o, p] = sum_value
    end
    return calc_zem_stiffness_tensor(C_rotated, Kinv)
end
