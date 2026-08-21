"""
    get_tensor(M, i)

$(extension_api_note())

Return column `i` of the storage field `M` as a `SMatrix{3,3}`, e.g. the deformation gradient
of point or bond `i` of a `PointTensor` or `BondTensor` field. The column is
read in column-major order, which is how [`update_tensor!`](@ref) writes it.

# Example

```julia
F = Peridynamics.get_tensor(storage.defgrad, i)
```

See also [`update_tensor!`](@ref), [`get_vector`](@ref).
"""
@inline function get_tensor(Mₙ::AbstractMatrix{T}, i::Int) where {T}
    tensor = SMatrix{3,3,T,9}(Mₙ[1, i], Mₙ[2, i], Mₙ[3, i], Mₙ[4, i], Mₙ[5, i], Mₙ[6, i],
                              Mₙ[7, i], Mₙ[8, i], Mₙ[9, i])
    return tensor
end

"""
    update_tensor!(M, i, T)

$(extension_api_note())

Write the `SMatrix{3,3}` `T` into column `i` of the storage field `M`, in column-major order.
The inverse of [`get_tensor`](@ref).

# Example

```julia
Peridynamics.update_tensor!(storage.defgrad, i, F)
```
"""
@inline function update_tensor!(Mₙ::AbstractMatrix{T}, i::Int,
                                Tₙ₊₁::StaticMatrix{3,3,T}) where {T}
    Mₙ[1, i] = Tₙ₊₁[1]
    Mₙ[2, i] = Tₙ₊₁[2]
    Mₙ[3, i] = Tₙ₊₁[3]
    Mₙ[4, i] = Tₙ₊₁[4]
    Mₙ[5, i] = Tₙ₊₁[5]
    Mₙ[6, i] = Tₙ₊₁[6]
    Mₙ[7, i] = Tₙ₊₁[7]
    Mₙ[8, i] = Tₙ₊₁[8]
    Mₙ[9, i] = Tₙ₊₁[9]
    return nothing
end

@inline function update_add_tensor!(Mₙ::AbstractMatrix{T}, i::Int,
                                    Tₙ₊₁::StaticMatrix{3,3,T}) where {T}
    Mₙ[1, i] += Tₙ₊₁[1]
    Mₙ[2, i] += Tₙ₊₁[2]
    Mₙ[3, i] += Tₙ₊₁[3]
    Mₙ[4, i] += Tₙ₊₁[4]
    Mₙ[5, i] += Tₙ₊₁[5]
    Mₙ[6, i] += Tₙ₊₁[6]
    Mₙ[7, i] += Tₙ₊₁[7]
    Mₙ[8, i] += Tₙ₊₁[8]
    Mₙ[9, i] += Tₙ₊₁[9]
    return nothing
end

@inline function zero_tensor!(Mₙ::AbstractMatrix{T}, i::Int) where {T}
    Mₙ[1, i] = zero(T)
    Mₙ[2, i] = zero(T)
    Mₙ[3, i] = zero(T)
    Mₙ[4, i] = zero(T)
    Mₙ[5, i] = zero(T)
    Mₙ[6, i] = zero(T)
    Mₙ[7, i] = zero(T)
    Mₙ[8, i] = zero(T)
    Mₙ[9, i] = zero(T)
    return nothing
end

"""
    get_sym_tensor(M, i)

$(extension_api_note())

Return column `i` of the storage field `M` as a symmetric `SMatrix{3,3}`, e.g. the plastic
strain of point or bond `i` of a `PointSymTensor` or `BondSymTensor` field.
The column holds the six independent components in Voigt order,
`(11, 22, 33, 23, 13, 12)`, which is how [`update_sym_tensor!`](@ref) writes it.

Note that a symmetric field has six rows, not nine, so [`get_tensor`](@ref) must not be used
on it.

# Example

```julia
εᵖ = Peridynamics.get_sym_tensor(state.bond_plastic_strain, idx)
```

See also [`update_sym_tensor!`](@ref), [`get_tensor`](@ref).
"""
@inline function get_sym_tensor(Mₙ::AbstractMatrix{T}, i::Int) where {T}
    tensor = SMatrix{3,3,T,9}(Mₙ[1, i], Mₙ[6, i], Mₙ[5, i],
                              Mₙ[6, i], Mₙ[2, i], Mₙ[4, i],
                              Mₙ[5, i], Mₙ[4, i], Mₙ[3, i])
    return tensor
end

"""
    update_sym_tensor!(M, i, T)

$(extension_api_note())

Write the symmetric `SMatrix{3,3}` `T` into column `i` of the storage field `M`, as the six
independent components in Voigt order, `(11, 22, 33, 23, 13, 12)`. The inverse of
[`get_sym_tensor`](@ref).

Only the upper triangle of `T` is read, and [`get_sym_tensor`](@ref) mirrors it back. A tensor
that is symmetric only up to round-off therefore comes back changed by that round-off, and a
tensor with a real skew part is *not* symmetrized — its lower triangle is silently discarded.

# Example

```julia
Peridynamics.update_sym_tensor!(state.bond_plastic_strain, idx, εᵖ + Δεᵖ)
```
"""
@inline function update_sym_tensor!(Mₙ::AbstractMatrix{T}, i::Int,
                                    Tₙ₊₁::StaticMatrix{3,3,T}) where {T}
    Mₙ[1, i] = Tₙ₊₁[1, 1]
    Mₙ[2, i] = Tₙ₊₁[2, 2]
    Mₙ[3, i] = Tₙ₊₁[3, 3]
    Mₙ[4, i] = Tₙ₊₁[2, 3]
    Mₙ[5, i] = Tₙ₊₁[1, 3]
    Mₙ[6, i] = Tₙ₊₁[1, 2]
    return nothing
end

"""
    get_vector(M, i)

$(extension_api_note())

Return column `i` of the storage field `M` as a `SVector{3}`, e.g. the position of point `i`
of a `PointVector` field.

# Example

```julia
u = Peridynamics.get_vector(storage.displacement, i)
```

See also [`update_vector!`](@ref), [`get_vector_diff`](@ref).
"""
@inline function get_vector(M::AbstractMatrix{T}, i::Int) where {T}
    return SVector{3,T}(M[1, i], M[2, i], M[3, i])
end

"""
    update_vector!(M, i, V)

$(extension_api_note())

Write the `SVector{3}` `V` into column `i` of the storage field `M`, overwriting what is
there. The inverse of [`get_vector`](@ref).

Note that a force density is accumulated over the bonds of a point, so it is written with
[`update_add_vector!`](@ref) and not with this function.
"""
@inline function update_vector!(Mₙ::AbstractMatrix{T}, i::Int,
                                Vₙ₊₁::StaticVector{3,T}) where {T}
    Mₙ[1, i] = Vₙ₊₁[1]
    Mₙ[2, i] = Vₙ₊₁[2]
    Mₙ[3, i] = Vₙ₊₁[3]
    return nothing
end

"""
    update_add_vector!(M, i, V)

$(extension_api_note())

Add the `SVector{3}` `V` to column `i` of the storage field `M`. This is how a force density
is accumulated inside `force_density_point!`, where every bond of a point contributes
a share.

# Example

```julia
Peridynamics.update_add_vector!(storage.b_int, i, b)
```
"""
@inline function update_add_vector!(Mₙ::AbstractMatrix{T}, i::Int,
                                    Vₙ₊₁::StaticVector{3,T}) where {T}
    Mₙ[1, i] += Vₙ₊₁[1]
    Mₙ[2, i] += Vₙ₊₁[2]
    Mₙ[3, i] += Vₙ₊₁[3]
    return nothing
end

"""
    get_vector_diff(M, i, j)

$(extension_api_note())

Return `column j - column i` of the storage field `M` as a `SVector{3}`, without building the
two columns first. This is the bond vector of a bond from point `i` to point `j`:

```julia
ΔXij = Peridynamics.get_vector_diff(system.position, i, j)  # initial bond vector
Δxij = Peridynamics.get_vector_diff(storage.position, i, j) # current bond vector
```
"""
@inline function get_vector_diff(M::AbstractMatrix{T}, i::Int, j::Int) where {T}
    return SVector{3,T}(M[1, j] - M[1, i], M[2, j] - M[2, i], M[3, j] - M[3, i])
end

"""
    invreg(M::StaticMatrix{N,N,T}, λ::Real, β::Real) where {N,T}

$(internal_api_warning())

Computes the regularized pseudo-inverse of a square static matrix `M` using a combination of
Tikhonov regularization in the SVD domain and truncated singular value regularization.

# Arguments

- `M::StaticMatrix{N,N,T}`: The square `N×N` static matrix with element type `T` to be
    inverted.
- `λ::Real`: Relative Tikhonov regularization parameter (dimensionless, non-negative).
    Controls the smoothing strength applied as
    ``\\lambda_{\\text{eff}} = \\lambda \\sigma_{\\max}``, where ``\\sigma_{\\max}`` is the
    largest singular value of `M`.
- `β::Real`: Relative SVD truncation parameter (dimensionless, non-negative). Defines the
    cutoff threshold as ``\\beta_{\\text{eff}} = \\beta \\sigma_{\\max}`` for excluding
    small singular values.

# Returns

- `Minv::StaticMatrix{N,N,T}`: The regularized pseudo-inverse of the input matrix `M`.

# Regularization Techniques

The function applies two complementary regularization strategies:

1. **SVD-based Tikhonov Regularization**: For each singular value ``\\sigma_i``, the inverse
    is computed as ``\\sigma_i/(\\sigma_i^2 + \\lambda_{\\text{eff}}^2)``, which smoothly
    dampens the contribution of small singular values without completely removing them.

2. **Truncated SVD**: Singular values below the threshold ``\\beta_{\\text{eff}}`` are
    completely excluded by setting their contribution to zero, preventing numerical
    instability from near-zero singular values.

!!! note "Scale-invariant regularization"
    Both ``\\lambda`` and ``\\beta`` are internally scaled by the largest singular value
    ``\\sigma_{\\max}``, making them **relative** regularization strengths independent of
    the matrix scale. This makes parameter selection more robust and transferable across
    different problems with varying magnitudes.

# Parameter Selection Guidelines

- **``λ`` (Tikhonov parameter)**:
    - Well-conditioned matrices: ``\\lambda = 0`` (no Tikhonov regularization, recommended
        default)
    - Mild regularization: ``\\lambda \\in [0, 10^{-12}]`` (scale-invariant gentle
        smoothing)
    - Moderate regularization: ``\\lambda \\in [10^{-12}, 10^{-4}]`` (for moderately
        ill-conditioned problems)
    - Note: Values ``\\lambda > 10^{-4}`` may introduce noticeable bias in the solution

- **``β`` (truncation parameter)**: Primary regularization mechanism, less sensitive than
    ``\\lambda``.
    - Well-conditioned matrices: ``\\beta \\in [\\sqrt{\\epsilon}, 10^{-6}]`` (remove
        numerical noise, recommended default)
    - Moderately ill-conditioned: ``\\beta \\in [10^{-6}, 10^{-4}]`` (moderate truncation)
    - Severely ill-conditioned: ``\\beta \\in [10^{-4}, 10^{-2}]`` (aggressive truncation)

See also the single-parameter method [`invreg(M, ε)`](@ref invreg), which regularizes each
singular value on its own and is therefore free of bias for the singular values that do not
need regularizing.
"""
function invreg(M::StaticMatrix{N,N,T}, λ::Real, β::Real) where {N,T}
    U, S, V = svd(M)
    λ_eff = λ * S[1] # the first singular value is the maximum
    β_eff = β * S[1]
    Sinvreg = SVector{N,T}((s > β_eff ? s/(s * s + λ_eff * λ_eff) : zero(T)) for s in S)
    Sinv = Diagonal{T,SVector{N,T}}(Sinvreg)
    return V * Sinv * U'
end

"""
    invreg(M::StaticMatrix{N,N,T}, ε::Real) where {N,T}

$(internal_api_warning())

Computes the regularized inverse of a square static matrix `M` with an adaptive,
per-singular-value Tikhonov damping that is only active where it is needed.

# Arguments

- `M::StaticMatrix{N,N,T}`: The square `N×N` static matrix with element type `T` to be
    inverted.
- `ε::Real`: Relative singular value floor (dimensionless, non-negative). The smallest
    singular value that is still trusted, expressed as a fraction of the largest one:
    ``\\sigma_f = \\varepsilon \\sigma_{\\max}``.

# Returns

- `Minv::StaticMatrix{N,N,T}`: The regularized inverse of the input matrix `M`.

# Regularization

Every singular value gets its own damping parameter

```math
\\lambda_i = \\max(0, \\sigma_f - \\sigma_i) , \\qquad
\\sigma_i^{\\mathrm{inv}} = \\frac{\\sigma_i}{\\sigma_i^2 + \\lambda_i^2}
```

so a singular value above the floor is damped with ``\\lambda_i = 0``, which is to say it is
inverted exactly, while one below the floor is damped smoothly towards zero:

| regime | ``\\lambda_i`` | ``\\sigma_i^{\\mathrm{inv}}`` |
|:---|:---|:---|
| ``\\sigma_i \\geq \\sigma_f`` (resolved) | ``0`` | ``1/\\sigma_i``, exactly |
| ``\\sigma_i = \\sigma_f / 2`` | ``\\sigma_f/2`` | ``1/\\sigma_f`` |
| ``\\sigma_i \\to 0`` (rank deficient) | ``\\sigma_f`` | ``\\to 0`` |

This bounds the inverse: the damped branch peaks at ``\\sigma_i = \\sigma_f/\\sqrt{2}``, so

```math
\\|M^{-1}\\|_2 \\leq \\frac{1}{(2\\sqrt{2}-2)\\,\\sigma_f}
              = \\frac{1}{(2\\sqrt{2}-2)\\,\\varepsilon\\,\\sigma_{\\max}}
              \\approx \\frac{1.21}{\\varepsilon\\,\\sigma_{\\max}}
```

no matter how rank deficient `M` becomes, and it stays bounded without a jump: the mapping is
``C^1`` across ``\\sigma_i = \\sigma_f``, because the
derivative of ``\\sigma/(\\sigma^2 + (\\sigma_f - \\sigma)^2)`` at ``\\sigma_f^-`` is
``-1/\\sigma_f^2``, which is the derivative of ``1/\\sigma`` at ``\\sigma_f^+``. A hard
truncation, as used by the two-parameter method above, is only ``C^0``: as a singular value
crosses the cutoff — which is what happens when bonds drop out of a family — the inverse jumps,
and a time integrator sees that jump as a force impulse.

!!! note "Relation to the two-parameter method"
    This method subsumes both ``\\lambda`` and ``\\beta`` of
    [`invreg(M, λ, β)`](@ref invreg) into a single parameter. A global ``\\lambda`` biases
    every singular value, and a global ``\\beta`` has to be small enough not to discard
    resolved ones, so no single pair can be both unbiased in the bulk and robust in the
    rank-deficient limit. Here the bias is exactly zero above the floor, which is what makes a
    floor as large as ``\\varepsilon = 10^{-3}`` affordable.

# Parameter Selection Guidelines

- ``\\varepsilon = 10^{-3}`` (recommended default): a 3D moment matrix of an intact point has
    ``\\sigma_{\\min}/\\sigma_{\\max} \\gtrsim 0.05``, so the regularization stays inactive in
    the bulk and only starts to act once enough bonds have failed to depress the conditioning
    by roughly 50x.
- ``\\varepsilon = 10^{-2}``: more conservative, and also acts on the intact matrices of thin
    structures with few points across the thickness.
- ``\\varepsilon = 0``: no regularization, the exact pseudo-inverse.
"""
function invreg(M::StaticMatrix{N,N,T}, ε::Real) where {N,T}
    U, S, V = svd(M)
    σ_f = ε * S[1] # the first singular value is the maximum
    Sinvreg = SVector{N,T}(begin
                               λ_i = max(zero(T), σ_f - s)
                               denom = s * s + λ_i * λ_i
                               # denom is only zero if s and σ_f both are, i.e. if M is zero
                               iszero(denom) ? zero(T) : s / denom
                           end for s in S)
    Sinv = Diagonal{T,SVector{N,T}}(Sinvreg)
    return V * Sinv * U'
end
