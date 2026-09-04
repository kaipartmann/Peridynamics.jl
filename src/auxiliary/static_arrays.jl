# The Voigt order used by `get_sym_tensor`/`update_sym_tensor!`: the diagonal first, then the
# off-diagonal pairs in the order the 3D convention `(11,22,33,23,13,12)` implies, i.e. the
# pair missing index `m` for `m` from 1 to `N` (only well-defined for `N` in (2, 3)).
function voigt_pairs(N::Int)
    if N == 2
        return ((1, 1), (2, 2), (1, 2))
    elseif N == 3
        return ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))
    else
        error("`get_sym_tensor`/`update_sym_tensor!` support N = 2 or N = 3, got N = $(N)!")
    end
end

# Voigt row index for every entry of the full `N × N` tensor, i.e. the inverse mapping of
# `voigt_pairs`, used to mirror the symmetric tensor back out in `get_sym_tensor`.
function voigt_rows(N::Int)
    pairs = voigt_pairs(N)
    return ntuple(r -> ntuple(c -> findfirst(==(minmax(r, c)), pairs), N), N)
end

"""
    get_tensor(M, i, ::Val{N})

$(extension_api_note())

Return column `i` of the storage field `M` as a `SMatrix{N,N}`, e.g. the deformation
gradient of point or bond `i` of a `PointTensor` or `BondTensor` field. The column is read
in column-major order, which is how [`update_tensor!`](@ref) writes it.

The number of spatial dimensions `N` is named by the `Val` of the last argument. It is
produced by [`dims`](@ref) from whatever is in scope, a system inside a force density kernel
and a storage inside a constitutive model hook.

# Example

```julia
F = Peridynamics.get_tensor(storage.defgrad, i, Peridynamics.dims(system))
```

See also [`update_tensor!`](@ref), [`get_vector`](@ref), [`dims`](@ref).
"""
@generated function get_tensor(M::AbstractMatrix{T}, i::Int, ::Val{N}) where {T,N}
    entries = [:(M[$d, i]) for d in 1:(N * N)]
    return quote
        @inline
        SMatrix{N,N,T,N * N}($(entries...))
    end
end

"""
    update_tensor!(M, i, A, ::Val{N})

$(extension_api_note())

Write the `SMatrix{N,N}` `A` into column `i` of the storage field `M`, in column-major
order. The inverse of [`get_tensor`](@ref). A value that is not a `StaticMatrix{N,N}` of the
`N` of the last argument is a `MethodError`, so a write can never disagree with the
dimension it names.

# Example

```julia
Peridynamics.update_tensor!(storage.defgrad, i, F, Peridynamics.dims(system))
```

See also [`get_tensor`](@ref), [`dims`](@ref).
"""
@generated function update_tensor!(Mₙ::AbstractMatrix{T}, i::Int,
                                   Aₙ₊₁::StaticMatrix{N,N,T}, ::Val{N}) where {T,N}
    stores = [:(Mₙ[$d, i] = Aₙ₊₁[$d]) for d in 1:(N * N)]
    return quote
        @inline
        $(stores...)
        return nothing
    end
end

@generated function update_add_tensor!(Mₙ::AbstractMatrix{T}, i::Int,
                                       Aₙ₊₁::StaticMatrix{N,N,T}, ::Val{N}) where {T,N}
    stores = [:(Mₙ[$d, i] += Aₙ₊₁[$d]) for d in 1:(N * N)]
    return quote
        @inline
        $(stores...)
        return nothing
    end
end

@generated function zero_tensor!(Mₙ::AbstractMatrix{T}, i::Int, ::Val{N}) where {T,N}
    stores = [:(Mₙ[$d, i] = zero(T)) for d in 1:(N * N)]
    return quote
        @inline
        $(stores...)
        return nothing
    end
end

"""
    get_sym_tensor(M, i, ::Val{N})

$(extension_api_note())

Return column `i` of the storage field `M` as a symmetric `SMatrix{N,N}`, e.g. the plastic
strain of point or bond `i` of a `PointSymTensor` or `BondSymTensor` field. The column holds
the independent components in Voigt order, the diagonal first and then the off-diagonals in
the order the 3D convention implies, `(11, 22, 33, 23, 13, 12)` for `N = 3` and
`(11, 22, 12)` for `N = 2`, which is how [`update_sym_tensor!`](@ref) writes it.

Note that a symmetric field has as many rows as there are independent components, not `N*N`,
so [`get_tensor`](@ref) must not be used on it.

The number of spatial dimensions `N` is named by the `Val` of the last argument, which
[`dims`](@ref) produces from whatever is in scope.

# Example

```julia
εᵖ = Peridynamics.get_sym_tensor(state.bond_plastic_strain, idx, Peridynamics.dims(storage))
```

See also [`update_sym_tensor!`](@ref), [`get_tensor`](@ref), [`dims`](@ref).
"""
@generated function get_sym_tensor(Mₙ::AbstractMatrix{T}, i::Int, ::Val{N}) where {T,N}
    rows = voigt_rows(N)
    entries = [:(Mₙ[$(rows[r][c]), i]) for c in 1:N for r in 1:N]
    return quote
        @inline
        SMatrix{N,N,T,N * N}($(entries...))
    end
end

"""
    update_sym_tensor!(M, i, A, ::Val{N})

$(extension_api_note())

Write the symmetric `SMatrix{N,N}` `A` into column `i` of the storage field `M`, as the
independent components in Voigt order, see [`get_sym_tensor`](@ref) for the order. The
inverse of [`get_sym_tensor`](@ref).

Only the upper triangle of `A` is read, and [`get_sym_tensor`](@ref) mirrors it back. A tensor
that is symmetric only up to round-off therefore comes back changed by that round-off, and a
tensor with a real skew part is not symmetrized. Its lower triangle is silently discarded.

# Example

```julia
Peridynamics.update_sym_tensor!(state.bond_plastic_strain, idx, εᵖ + Δεᵖ,
                                Peridynamics.dims(storage))
```

See also [`get_sym_tensor`](@ref), [`dims`](@ref).
"""
@generated function update_sym_tensor!(Mₙ::AbstractMatrix{T}, i::Int,
                                       Aₙ₊₁::StaticMatrix{N,N,T}, ::Val{N}) where {T,N}
    pairs = voigt_pairs(N)
    stores = [:(Mₙ[$k, i] = Aₙ₊₁[$r, $c]) for (k, (r, c)) in enumerate(pairs)]
    return quote
        @inline
        $(stores...)
        return nothing
    end
end

"""
    get_vector(M, i, ::Val{N})

$(extension_api_note())

Return column `i` of the storage field `M` as a `SVector{N}`, e.g. the position of point `i`
of a `PointVector` field.

The number of spatial dimensions `N` is named by the `Val` of the last argument, which
[`dims`](@ref) produces from whatever is in scope.

# Example

```julia
u = Peridynamics.get_vector(storage.displacement, i, Peridynamics.dims(system))
```

See also [`update_vector!`](@ref), [`get_vector_diff`](@ref), [`dims`](@ref).
"""
@generated function get_vector(M::AbstractMatrix{T}, i::Int, ::Val{N}) where {T,N}
    entries = [:(M[$d, i]) for d in 1:N]
    return quote
        @inline
        SVector{N,T}($(entries...))
    end
end

"""
    update_vector!(M, i, V, ::Val{N})

$(extension_api_note())

Write the `SVector{N}` `V` into column `i` of the storage field `M`, overwriting what is
there. The inverse of [`get_vector`](@ref). A value that is not a `StaticVector{N}` of the
`N` of the last argument is a `MethodError`.

Note that a force density is accumulated over the bonds of a point, so it is written with
[`update_add_vector!`](@ref) and not with this function.

See also [`get_vector`](@ref), [`dims`](@ref).
"""
@generated function update_vector!(Mₙ::AbstractMatrix{T}, i::Int,
                                   Vₙ₊₁::StaticVector{N,T}, ::Val{N}) where {N,T}
    stores = [:(Mₙ[$d, i] = Vₙ₊₁[$d]) for d in 1:N]
    return quote
        @inline
        $(stores...)
        return nothing
    end
end

"""
    update_add_vector!(M, i, V, ::Val{N})

$(extension_api_note())

Add the `SVector{N}` `V` to column `i` of the storage field `M`. This is how a force density
is accumulated inside `force_density_point!`, where every bond of a point contributes a
share.

# Example

```julia
Peridynamics.update_add_vector!(storage.b_int, i, b, Peridynamics.dims(system))
```

See also [`update_vector!`](@ref), [`dims`](@ref).
"""
@generated function update_add_vector!(Mₙ::AbstractMatrix{T}, i::Int,
                                       Vₙ₊₁::StaticVector{N,T}, ::Val{N}) where {N,T}
    stores = [:(Mₙ[$d, i] += Vₙ₊₁[$d]) for d in 1:N]
    return quote
        @inline
        $(stores...)
        return nothing
    end
end

"""
    get_vector_diff(M, i, j, ::Val{N})

$(extension_api_note())

Return `column j - column i` of the storage field `M` as a `SVector{N}`, without building
the two columns first. This is the bond vector of a bond from point `i` to point `j`:

```julia
# initial bond vector
ΔXij = Peridynamics.get_vector_diff(system.position, i, j, Peridynamics.dims(system))
# current bond vector
Δxij = Peridynamics.get_vector_diff(storage.position, i, j, Peridynamics.dims(system))
```

See also [`get_vector`](@ref), [`dims`](@ref).
"""
@generated function get_vector_diff(M::AbstractMatrix{T}, i::Int, j::Int,
                                    ::Val{N}) where {T,N}
    entries = [:(M[$d, j] - M[$d, i]) for d in 1:N]
    return quote
        @inline
        SVector{N,T}($(entries...))
    end
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
crosses the cutoff, which is what happens when bonds drop out of a family, the inverse jumps,
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
