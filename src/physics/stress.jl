function von_mises_stress(σ)
    σx, σy, σz = σ[1], σ[5], σ[9]
    τxy, τxz, τyz = σ[4], σ[7], σ[8]
    a = σx * σx + σy * σy + σz * σz
    b = - σx * σy - σx * σz - σy * σz
    c = 3 * (τxy * τxy + τxz * τxz + τyz * τyz)
    d = a + b + c
    d < 0 && return zero(eltype(σ))
    σvm = √d
    return σvm
end

function cauchy_stress(P, F)
    J = det(F)
    σ = 1/J .* P * F'
    return σ
end

"""
    cauchy_stress_safe(P, F)

Compute the Cauchy stress tensor from the first Piola-Kirchhoff stress and deformation
gradient with NaN protection. Returns zero stress if the Jacobian is too small or NaN.
"""
function cauchy_stress_safe(P, F)
    J = det(F)
    if J < 1e-30 || isnan(J)
        return zero(SMatrix{3,3,Float64,9})
    end
    σ = 1/J .* P * F'
    return σ
end

function first_piola_kirchhoff(σ, F)
    J = det(F)
    P = J * σ * inv(F)'
    return P
end

function init_stress_rotation!(storage::AbstractStorage, F, Ḟ, Δt, i)
    # inverse of the deformation gradient
    F⁻¹ = inv(F)
    if containsnan(F⁻¹)
        OTensor = zero(SMatrix{3,3,Float64,9})
        update_tensor!(storage.rotation, i, OTensor)
        update_tensor!(storage.left_stretch, i, OTensor)
        return OTensor
    end

    # Eulerian velocity gradient [FT87, eq. (3)]
    L = Ḟ * F⁻¹

    # rate-of-deformation tensor D
    D = 0.5 .* (L + L')

    # spin tensor W
    W = 0.5 .* (L - L')

    # left stretch V
    V = get_tensor(storage.left_stretch, i)

    # vector z [FT87, eq. (13)]
    z_x = - V[1,3] * D[2,1] - V[2,3] * D[2,2] -
            V[3,3] * D[2,3] + V[1,2] * D[3,1] +
            V[2,2] * D[3,2] + V[3,2] * D[3,3]
    z_y = V[1,3] * D[1,1] + V[2,3] * D[1,2] +
          V[3,3] * D[1,3] - V[1,1] * D[3,1] -
          V[2,1] * D[3,2] - V[3,1] * D[3,3]
    z_z = - V[1,2] * D[1,1] - V[2,2] * D[1,2] -
            V[3,2] * D[1,3] + V[1,1] * D[2,1] +
            V[2,1] * D[2,2] + V[3,1] * D[2,3]
    z = SVector{3}(z_x, z_y, z_z)

    # w = -1/2 * \epsilon_{ijk} * W_{jk}  [FT87, eq. (11)]
    w = 0.5 .* SVector{3}(W[3,2] - W[2,3], W[1,3] - W[3,1], W[2,1] - W[1,2])

    # ω = w + (I * tr(V) - V)^(-1) * z [FT87, eq. (12)]
    ω = w + inv(I * tr(V) - V) * z

    # Ω [FT87, eq. (10)]
    Ωt1, Ωt4, Ωt7 = 0.0, -ω[3], ω[2]
    Ωt2, Ωt5, Ωt8 = ω[3], 0.0, -ω[1]
    Ωt3, Ωt6, Ωt9 = -ω[2], ω[1], 0.0
    Ωtens = SMatrix{3,3}(Ωt1, Ωt2, Ωt3, Ωt4, Ωt5, Ωt6, Ωt7, Ωt8, Ωt9)
    Ω² = dot(ω, ω)
    Ω = sqrt(Ω²)

    # compute Q with [FT87, eq. (44)]
    if Ω² > 1e-30 # avoid a potential divide-by-zero
        fac1 = sin(Δt * Ω) / Ω
        fac2 = -(1.0 - cos(Δt * Ω)) / Ω²
        Ωtens² = Ωtens * Ωtens
        Q = I + fac1 .* Ωtens + fac2 .* Ωtens²
    else
        Q = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    end

    # compute Rotation of new step [FT87, eq. (36)]
    R = get_tensor(storage.rotation, i)
    Rₙ₊₁ = Q * R

    # compute step 4 of [FT87]
    V̇ = L * V - V * Ωtens

    # compute step 5 of [FT87]
    Vₙ₊₁ = V + Δt * V̇

    # compute the unrotated rate of deformation
    Dᵣ = Rₙ₊₁' * D * Rₙ₊₁

    # update rotation and left stretch
    update_tensor!(storage.rotation, i, Rₙ₊₁)
    update_tensor!(storage.left_stretch, i, Vₙ₊₁)
    return Dᵣ
end

function rotate_stress(storage::AbstractStorage, σ, i)
    R = get_tensor(storage.rotation, i)
    T = R * σ * R'
    return T
end

#-------------------------------------------------------------------------------------------
# Phase-field style stress splitting functions
#-------------------------------------------------------------------------------------------

"""
    stress_split_spectral(P::SMatrix{3,3}, F::SMatrix{3,3})

Compute the spectral decomposition of the first Piola-Kirchhoff stress tensor into
tensile (positive) and compressive (negative) parts. This follows the phase-field
fracture approach by Miehe et al. (2010).

Returns a tuple `(P_plus, P_minus, success)` where:
- `P_plus`: The tensile part of the stress (associated with positive principal stresses)
- `P_minus`: The compressive part of the stress (associated with negative principal stresses)
- `success`: Boolean indicating whether the decomposition was successful

The split is performed in Cauchy stress space using spectral decomposition, then
converted back to first Piola-Kirchhoff stress.

# Phase-field energy splitting background

In phase-field fracture, the strain energy is split as:
```math
\\Psi(\\epsilon, d) = g(d) \\Psi^+(\\epsilon) + \\Psi^-(\\epsilon)
```
where `g(d) = (1-d)²` is the degradation function and only the tensile part `Ψ⁺`
drives fracture. This prevents unrealistic crack growth under compression.

# Reference
Miehe, C., Hofacker, M., & Welschinger, F. (2010). A phase field model for rate-independent
crack propagation. *Computer Methods in Applied Mechanics and Engineering*, 199(45-48), 2765-2778.
"""
function stress_split_spectral(P::SMatrix{3,3,T,9}, F::SMatrix{3,3,T,9}) where {T}
    zero_P = zero(SMatrix{3,3,T,9})

    # Check deformation gradient validity
    J = det(F)
    if J < 1e-10 || isnan(J) || isinf(J)
        return zero_P, zero_P, false
    end

    # Check if stress is already very small or has NaN
    P_norm = norm(P)
    if isnan(P_norm) || isinf(P_norm)
        return zero_P, zero_P, false
    end
    if P_norm < 1e-30
        return zero_P, zero_P, true  # zero stress splits into zeros
    end

    # Convert to Cauchy stress for spectral decomposition
    σ = (1/J) .* P * F'

    # Check for NaN/Inf in Cauchy stress
    if any(isnan, σ) || any(isinf, σ)
        return zero_P, zero_P, false
    end

    # Spectral decomposition of Cauchy stress
    σ_sym = Symmetric(σ)
    eig_result = eigen(σ_sym)
    λ = eig_result.values
    V = eig_result.vectors

    # Check eigendecomposition validity
    if any(isnan, λ) || any(isinf, λ) || any(isnan, V) || any(isinf, V)
        return zero_P, zero_P, false
    end

    # Split into positive and negative parts using Macaulay brackets
    σ_plus = zero(SMatrix{3,3,T,9})
    σ_minus = zero(SMatrix{3,3,T,9})

    for i in 1:3
        ni = SVector{3,T}(V[1,i], V[2,i], V[3,i])
        proj = ni * ni'  # outer product (projection tensor)
        if λ[i] > 0
            σ_plus += λ[i] * proj
        else
            σ_minus += λ[i] * proj
        end
    end

    # Convert back to first Piola-Kirchhoff stress using F^{-T}
    # Use pseudo-inverse approach for numerical stability
    F_inv = try
        inv(F)
    catch
        return zero_P, zero_P, false
    end

    if any(isnan, F_inv) || any(isinf, F_inv)
        return zero_P, zero_P, false
    end

    Finv_T = F_inv'
    P_plus = J * σ_plus * Finv_T
    P_minus = J * σ_minus * Finv_T

    # Final validity check
    if any(isnan, P_plus) || any(isinf, P_plus) || any(isnan, P_minus) || any(isinf, P_minus)
        return zero_P, zero_P, false
    end

    return P_plus, P_minus, true
end

"""
    degraded_stress(P::SMatrix{3,3}, F::SMatrix{3,3}, d::Real)

Apply phase-field style degradation to the stress tensor with tensile/compressive
splitting. The degradation function `g(d) = (1-d)²` is applied only to the tensile
part of the stress, while the compressive part remains undegraded.

```math
\\boldsymbol{P}_{\\text{eff}} = g(d) \\boldsymbol{P}^+ + \\boldsymbol{P}^-
```

This prevents unrealistic behavior under compression and is consistent with the
physics of brittle fracture where cracks cannot propagate under pure compression.

If the spectral decomposition fails (due to ill-conditioned deformation gradient),
falls back to standard degradation: `g(d) * P`.

A minimum stiffness floor is enforced: `g(d) = max(ε, (1-d)²)` to prevent complete
loss of stiffness and associated numerical instabilities.

# Arguments
- `P`: First Piola-Kirchhoff stress tensor
- `F`: Deformation gradient
- `d`: Damage variable (0 = intact, 1 = fully broken)

# Returns
The effective (degraded) first Piola-Kirchhoff stress tensor.
"""
function degraded_stress(P::SMatrix{3,3,T,9}, F::SMatrix{3,3,T,9}, d::Real) where {T}
    # Minimum stiffness floor to prevent complete loss of resistance
    # This improves numerical stability near full damage
    g_d_min = 1e-6
    g_d = max(g_d_min, (1.0 - d) * (1.0 - d))  # quadratic degradation with floor

    P_plus, P_minus, success = stress_split_spectral(P, F)

    if !success
        # Fallback to standard (non-split) degradation
        return g_d * P
    end

    return g_d * P_plus + P_minus
end
