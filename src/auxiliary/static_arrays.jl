
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
    V = SVector{3}(M[1, i], M[2, i], M[3, i])
    return V
end

@inline function update_vector!(Mₙ::AbstractMatrix, i::Int, Vₙ₊₁::SVector{3})
    Mₙ[1, i] = Vₙ₊₁[1]
    Mₙ[2, i] = Vₙ₊₁[2]
    Mₙ[3, i] = Vₙ₊₁[3]
    return nothing
end

@inline function get_vector_diff(M::AbstractMatrix, i::Int, j::Int)
    V = SVector{3}(M[1, j] - M[1, i], M[2, j] - M[2, i], M[3, j] - M[3, i])
    return V
end
