function von_mises_stress(σ)
    σx, σy, σz = σ[1], σ[5], σ[9]
    τxy, τxz, τyz = σ[4], σ[7], σ[8]
    a = σx * σx + σy * σy + σz * σz
    b = - σx * σy - σx * σz - σy * σz
    c = 3 * (τxy * τxy + τxz * τxz + τyz * τyz)
    σvm = √(a + b + c)
    return σvm
end

function cauchy_stress(P, F)
    J = det(F)
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
    Ωtens = SMatrix{3,3}(0.0, -ω[3], ω[2], ω[3], 0.0, -ω[1], -ω[2], ω[1], 0.0)
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
