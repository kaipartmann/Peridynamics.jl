"""
    sym_eigvals(A)

$(extension_api_note())

Return the three eigenvalues of the symmetric `SMatrix{3,3}` `A` in descending order, as an
`SVector{3}`. They are computed in closed form (Smith, 1961) rather than by an
eigendecomposition, so this allocates nothing and needs no eigenvectors.

Only the upper triangle of `A` is read.

See also [`hencky_and_invstretch`](@ref).
"""
@inline function sym_eigvals(A::StaticMatrix{3,3,T}) where {T}
    p1 = A[1, 2]^2 + A[1, 3]^2 + A[2, 3]^2
    if p1 < eps(T)
        # a diagonal tensor: the diagonal is the spectrum, but it is not necessarily sorted
        return sort_desc(A[1, 1], A[2, 2], A[3, 3])
    end
    q = tr(A) / 3
    p2 = (A[1, 1] - q)^2 + (A[2, 2] - q)^2 + (A[3, 3] - q)^2 + 2 * p1
    p = sqrt(p2 / 6)
    B = (A - q * I) / p
    r = clamp(det(B) / 2, -one(T), one(T))
    φ = acos(r) / 3
    λ1 = q + 2 * p * cos(φ)
    λ3 = q + 2 * p * cos(φ + 2 * T(π) / 3)
    λ2 = 3 * q - λ1 - λ3
    return SVector{3,T}(λ1, λ2, λ3)
end

@inline function sort_desc(a::T, b::T, c::T) where {T}
    a < b && ((a, b) = (b, a))
    b < c && ((b, c) = (c, b))
    a < b && ((a, b) = (b, a))
    return SVector{3,T}(a, b, c)
end

"""
    hencky_and_invstretch(C)

$(extension_api_note())

Return the material logarithmic (Hencky) strain ``\\boldsymbol{\\varepsilon} = \\ln
\\boldsymbol{U} = \\frac{1}{2} \\ln \\boldsymbol{C}`` and the inverse right stretch
``\\boldsymbol{U}^{-1} = \\boldsymbol{C}^{-1/2}`` of the right Cauchy-Green tensor
``\\boldsymbol{C} = \\boldsymbol{F}^T \\boldsymbol{F}``, as a tuple of two `SMatrix{3,3}`.

Both are isotropic functions of `C` and are evaluated in closed form, from the analytic
eigenvalues of [`sym_eigvals`](@ref) and the spectral projectors built from them by Lagrange
interpolation, so no eigenvectors are computed and nothing is allocated. Repeated eigenvalues
are handled by the reduced interpolation they call for.

This is the pair a finite-strain constitutive model needs to work in logarithmic strain
space: ``\\boldsymbol{\\varepsilon}`` and the rotated Kirchhoff stress
``\\hat{\\boldsymbol{\\tau}}`` are work conjugate, and the pull-back to the first
Piola-Kirchhoff stress is ``\\boldsymbol{P} = \\boldsymbol{F} \\, \\boldsymbol{U}^{-1}
\\hat{\\boldsymbol{\\tau}} \\, \\boldsymbol{U}^{-1}``.

# Example

```julia
ε, Uinv = Peridynamics.hencky_and_invstretch(F' * F)
τ = params.λ * tr(ε) * I + 2 * params.μ * ε
P = F * (Uinv * τ * Uinv)
```
"""
function hencky_and_invstretch(C::StaticMatrix{3,3,T}) where {T}
    c = sym_eigvals(C)
    c1, c2, c3 = max(c[1], eps(T)), max(c[2], eps(T)), max(c[3], eps(T))
    Id = SMatrix{3,3,T,9}(I)

    # eigenvalues closer than this are treated as repeated, because the Lagrange
    # interpolation below cancels catastrophically when they are not separated
    tol = sqrt(eps(T)) * max(c1, one(T))
    c12_equal = c1 - c2 < tol
    c23_equal = c2 - c3 < tol

    if c12_equal && c23_equal # spherical: one eigenvalue of multiplicity three
        return log(c1) / 2 * Id, 1 / sqrt(c1) * Id
    elseif c12_equal # c1 = c2 ≠ c3, so one projector is enough
        P3 = (C - c1 * Id) / (c3 - c1)
        return log(c3) / 2 * P3 + log(c1) / 2 * (Id - P3),
               1 / sqrt(c3) * P3 + 1 / sqrt(c1) * (Id - P3)
    elseif c23_equal # c1 ≠ c2 = c3
        P1 = (C - c2 * Id) / (c1 - c2)
        return log(c1) / 2 * P1 + log(c2) / 2 * (Id - P1),
               1 / sqrt(c1) * P1 + 1 / sqrt(c2) * (Id - P1)
    end

    # three distinct eigenvalues: the spectral projectors, built without eigenvectors
    P1 = ((C - c2 * Id) * (C - c3 * Id)) / ((c1 - c2) * (c1 - c3))
    P2 = ((C - c1 * Id) * (C - c3 * Id)) / ((c2 - c1) * (c2 - c3))
    P3 = Id - P1 - P2
    ε = log(c1) / 2 * P1 + log(c2) / 2 * P2 + log(c3) / 2 * P3
    Uinv = 1 / sqrt(c1) * P1 + 1 / sqrt(c2) * P2 + 1 / sqrt(c3) * P3
    return ε, Uinv
end

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
