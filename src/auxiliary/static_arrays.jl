
@inline function get_tensor(Mₙ::AbstractMatrix{T}, i::Int) where {T}
    tensor = SMatrix{3,3,T,9}(Mₙ[1, i], Mₙ[2, i], Mₙ[3, i], Mₙ[4, i], Mₙ[5, i], Mₙ[6, i],
                              Mₙ[7, i], Mₙ[8, i], Mₙ[9, i])
    return tensor
end

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

@inline function get_vector(M::AbstractMatrix{T}, i::Int) where {T}
    return SVector{3,T}(M[1, i], M[2, i], M[3, i])
end

@inline function update_vector!(Mₙ::AbstractMatrix{T}, i::Int,
                                Vₙ₊₁::StaticVector{3,T}) where {T}
    Mₙ[1, i] = Vₙ₊₁[1]
    Mₙ[2, i] = Vₙ₊₁[2]
    Mₙ[3, i] = Vₙ₊₁[3]
    return nothing
end

@inline function update_add_vector!(Mₙ::AbstractMatrix{T}, i::Int,
                                    Vₙ₊₁::StaticVector{3,T}) where {T}
    Mₙ[1, i] += Vₙ₊₁[1]
    Mₙ[2, i] += Vₙ₊₁[2]
    Mₙ[3, i] += Vₙ₊₁[3]
    return nothing
end

@inline function get_vector_diff(M::AbstractMatrix{T}, i::Int, j::Int) where {T}
    return SVector{3,T}(M[1, j] - M[1, i], M[2, j] - M[2, i], M[3, j] - M[3, i])
end

"""
    invreg(M::StaticMatrix{N,N,T}, λ::Real, β::Real) where {N,T}

Computes the regularized pseudo-inverse of a square static matrix `M` using a combination of
Tikhonov regularization in the SVD domain and truncated singular value regularization.

# Arguments

- `M::StaticMatrix{N,N,T}`: The square `N`×`N` static matrix with element type `T` to be
    inverted.
- `λ::Real`: Relative Tikhonov regularization parameter (dimensionless, non-negative).
    Controls the smoothing strength applied as λₑ = λ σₘₐₓ, where σₘₐₓ is the largest
    singular value of `M`.
- `β::Real`: Relative SVD truncation parameter (dimensionless, non-negative). Defines the
    cutoff threshold as βₑ = β σₘₐₓ for excluding small singular values.

# Returns

- `Minv::StaticMatrix{N,N,T}`: The regularized pseudo-inverse of the input matrix `M`.

# Regularization Techniques

The function applies two complementary regularization strategies:

1. **SVD-based Tikhonov Regularization**: For each singular value σᵢ, the inverse is
   computed as σᵢ/(σᵢ² + λₑ²), which smoothly dampens the contribution of small singular
   values without completely removing them.

2. **Truncated SVD**: Singular values below the threshold βₑ are completely excluded by
   setting their contribution to zero, preventing numerical instability from near-zero
   singular values.

!!! note "Scale-invariant regularization"
    Both λ and β are internally scaled by the largest singular value σₘₐₓ, making them
    **relative** regularization strengths independent of the matrix scale. This makes
    parameter selection more robust and transferable across different problems with
    varying magnitudes.

# Parameter Selection Guidelines

- **λ (Tikhonov parameter)**:
  - Well-conditioned matrices: λ = 0 (no Tikhonov regularization, recommended default)
  - Mild regularization: λ ∈ [0, 10⁻¹²] (scale-invariant gentle smoothing)
  - Moderate regularization: λ ∈ [10⁻¹², 10⁻⁴] (for moderately ill-conditioned problems)
  - Note: Values λ > 10⁻⁴ may introduce noticeable bias in the solution

- **β (truncation parameter)**: Primary regularization mechanism, less sensitive than λ.
  - Well-conditioned matrices: β ∈ [√ε, 10⁻⁶] (remove numerical noise, recommended default)
  - Moderately ill-conditioned: β ∈ [10⁻⁶, 10⁻⁴] (moderate truncation)
  - Severely ill-conditioned: β ∈ [10⁻⁴, 10⁻²] (aggressive truncation)
"""
function invreg(M::StaticMatrix{N,N,T}, λ::Real, β::Real) where {N,T}
    U, S, V = svd(M)
    λ_eff = λ * S[1] # the first singular value is the maximum
    β_eff = β * S[1]
    Sinvreg = SVector{N,T}((s > β_eff ? s/(s * s + λ_eff * λ_eff) : zero(T)) for s in S)
    Sinv = Diagonal{T,SVector{N,T}}(Sinvreg)
    return V * Sinv * U'
end
