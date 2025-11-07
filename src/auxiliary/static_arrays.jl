
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

Computes the pseudo-inverse of a square static matrix `M` using a combination of
Tikhonov regularization and SVD-based truncated singular value regularization.

# Regularization Techniques

1. **Tikhonov Regularization**: The matrix `M` is first regularized by adding
   `λ * tr(M) * I` to it, which helps improve conditioning by shifting eigenvalues
   away from zero.

2. **Truncated SVD**: After computing the singular value decomposition, singular
   values below `β * S[1]` (where `S[1]` is the largest singular value) are
   replaced by zero in the pseudo-inverse, preventing numerical instability from
   very small singular values.

# Parameter Selection Guidelines

- **Well-conditioned matrices**: Use small values (λ ∈ [0, 10⁻¹⁰], β ∈ [0, 10⁻⁸])
    to minimize regularization effects while maintaining numerical stability.
- **Moderately ill-conditioned matrices**: Use moderate values (λ ∈ [10⁻⁶, 10⁻³],
    β ∈ [10⁻⁸, 10⁻⁶]) to balance accuracy and stability.
- **Severely ill-conditioned matrices**: Use larger values (λ ∈ [10⁻³, 10⁻¹],
    β ∈ [10⁻⁶, 10⁻⁴]) to ensure numerical stability at the cost of some accuracy.

# Arguments

- `M::StaticMatrix{N,N,T}`: The square `N`×`N` static matrix with element type `T` to be
    inverted.
- `λ::Real`: Tikhonov regularization parameter. Controls the strength of the
    regularization applied to the matrix. Should be a non-negative real number.
    Larger values increase regularization strength.
- `β::Real`: SVD truncation parameter. Defines the threshold as a fraction of the
    largest singular value. Should be a positive real number between 0 and 1.
    Singular values below `β * S[1]` are set to zero in the pseudo-inverse.

# Returns

- `Minv::StaticMatrix{N,N,T}`: The regularized pseudo-inverse of the input matrix `M`.
"""
function invreg(M::StaticMatrix{N,N,T}, λ::Real, β::Real) where {N,T}
    # Tikhonov regularization
    Mreg = M + λ * tr(M) * I

    # Truncated SVD-based regularized inversion
    U, S, V = svd(Mreg)
    threshold = β * S[1] # the first singular value is the maximum
    Sinvreg = SVector{N,T}((s > threshold ? 1/s : zero(T)) for s in S)
    Sinv = Diagonal{T,SVector{N,T}}(Sinvreg)
    Minv = V * Sinv * U'

    return Minv
end
