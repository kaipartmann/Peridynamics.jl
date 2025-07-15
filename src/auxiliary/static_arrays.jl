
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
