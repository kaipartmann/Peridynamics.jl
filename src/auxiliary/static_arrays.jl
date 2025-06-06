
@inline function get_tensor(T::AbstractMatrix, i::Int)
    tensor = SMatrix{3,3}(T[1, i], T[2, i], T[3, i], T[4, i], T[5, i], T[6, i], T[7, i],
                          T[8, i], T[9, i])
    return tensor
end

@inline function update_tensor!(Tₙ::AbstractMatrix, i::Int, Tₙ₊₁::SMatrix{3,3})
    Tₙ[1, i] = Tₙ₊₁[1, 1]
    Tₙ[2, i] = Tₙ₊₁[1, 2]
    Tₙ[3, i] = Tₙ₊₁[1, 3]
    Tₙ[4, i] = Tₙ₊₁[2, 1]
    Tₙ[5, i] = Tₙ₊₁[2, 2]
    Tₙ[6, i] = Tₙ₊₁[2, 3]
    Tₙ[7, i] = Tₙ₊₁[3, 1]
    Tₙ[8, i] = Tₙ₊₁[3, 2]
    Tₙ[9, i] = Tₙ₊₁[3, 3]
    return nothing
end

@inline function get_vector(M::AbstractMatrix, i::Int)
    return SVector{3}(M[1, i], M[2, i], M[3, i])
end

@inline function update_vector!(Mₙ::AbstractMatrix, i::Int, Vₙ₊₁::SVector{3})
    Mₙ[1, i] = Vₙ₊₁[1]
    Mₙ[2, i] = Vₙ₊₁[2]
    Mₙ[3, i] = Vₙ₊₁[3]
    return nothing
end

@inline function update_add_vector!(Mₙ::AbstractMatrix, i::Int, Vₙ₊₁::SVector{3})
    Mₙ[1, i] += Vₙ₊₁[1]
    Mₙ[2, i] += Vₙ₊₁[2]
    Mₙ[3, i] += Vₙ₊₁[3]
    return nothing
end

@inline function get_vector_diff(M::AbstractMatrix, i::Int, j::Int)
    return SVector{3}(M[1, j] - M[1, i], M[2, j] - M[2, i], M[3, j] - M[3, i])
end

"""
    invreg(M::StaticMatrix{N,N,T}, threshold::Real=eps()) where {N,T}

Computes the pseudo-inverse of a square static matrix `M` using singular value
decomposition (SVD) with regularization. The regularization is applied to the
singular values, where singular values below the specified `threshold` are set
to zero in the pseudo-inverse. This helps to avoid numerical instability and
ill-conditioning in the inversion process.

# Arguments
- `M::StaticMatrix{N,N,T}`: The square `N`×`N` static matrix with element type `T` to be
    inverted.
- `threshold::Real=eps()`: The threshold value for regularization, should be a positive
    real number between 0 and 1. Singular values below this threshold are set to zero in the
    pseudo-inverse.
"""
function invreg(M::StaticMatrix{N,N,T}, threshold::Real=eps()) where {N,T}
    U, S, V = svd(M)
    Sinvreg = SVector{N,T}((s > threshold ? 1/s : 0) for s in S)
    Sinv = Diagonal{T,SVector{N,T}}(Sinvreg)
    Minv = V * Sinv * U'
    return Minv
end
