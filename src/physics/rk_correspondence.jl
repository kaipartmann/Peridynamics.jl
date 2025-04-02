"""
    RKCMaterial(; kernel, model, dmgmodel, maxdmg)

A material type used to assign the material of a [`Body`](@ref) with the local continuum
consistent (correspondence) formulation of non-ordinary state-based peridynamics.

# Keywords
- `kernel::Function`: Kernel function used for weighting the interactions between points. \\
    (default: `cubic_b_spline_kernel`) \\
    The following kernels can be used:
    - [`cubic_b_spline_kernel`](@ref)
- `model::AbstractConstitutiveModel`: Constitutive model defining the material behavior. \\
    (default: `LinearElastic()`) \\
    The following models can be used:
    - [`LinearElastic`](@ref)
    - [`NeoHooke`](@ref)
    - [`MooneyRivlin`](@ref)
    - [`SaintVenantKirchhoff`](@ref)
- `dmgmodel::AbstractDamageModel`: Damage model defining the damage behavior. \\
    (default: [`CriticalStretch`](@ref))
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. If this value is
    exceeded, all bonds of that point are broken because the deformation gradient would then
    possibly contain `NaN` values. \\
    (default: `0.85`)

!!! note "Stability of fracture simulations"
    This formulation is known to be not suitable for fracture simulations without
    stabilization of the zero-energy modes. Therefore be careful when doing fracture
    simulations and try out different parameters for `maxdmg` and `zem`.

# Examples

```julia-repl
julia> mat = RKCMaterial()
RKCMaterial{LinearElastic, typeof(linear_kernel), CriticalStretch}(maxdmg=0.85)
```

---

```julia
RKCMaterial{CM,K,DM}
```

Material type for the naturally stabilized local continuum consistent (correspondence) formulation of
non-ordinary state-based peridynamics.

# Type Parameters
- `CM`: A constitutive model type. See the constructor docs for more informations.
- `K`: A kernel function type. See the constructor docs for more informations.
- `DM`: A damage model type. See the constructor docs for more informations.

# Fields
- `kernel::Function`: Kernel function used for weighting the interactions between points.
    See the constructor docs for more informations.
- `model::AbstractConstitutiveModel`: Constitutive model defining the material behavior. See
    the constructor docs for more informations.
- `dmgmodel::AbstractDamageModel`: Damage model defining the damage behavior. See the
    constructor docs for more informations.
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. See the
    constructor docs for more informations.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `RKCMaterial`, then the following
parameters are allowed:
Material parameters:
- `horizon::Float64`: Radius of point interactions
- `rho::Float64`: Density
Elastic parameters:
- `E::Float64`: Young's modulus
- `nu::Float64`: Poisson's ratio
- `G::Float64`: Shear modulus
- `K::Float64`: Bulk modulus
- `lambda::Float64`: 1st Lamé parameter
- `mu::Float64`: 2nd Lamé parameter
Fracture parameters:
- `Gc::Float64`: Critical energy release rate
- `epsilon_c::Float64`: Critical strain

!!! note "Elastic parameters"
    Note that exactly two elastic parameters are required to specify a material.
    Please choose two out of the six allowed elastic parameters.

# Allowed export fields
When specifying the `fields` keyword of [`Job`](@ref) for a [`Body`](@ref) with
`RKCMaterial`, the following fields are allowed:
- `position::Matrix{Float64}`: Position of each point
- `displacement::Matrix{Float64}`: Displacement of each point
- `velocity::Matrix{Float64}`: Velocity of each point
- `velocity_half::Matrix{Float64}`: Velocity parameter for Verlet time solver
- `acceleration::Matrix{Float64}`: Acceleration of each point
- `b_int::Matrix{Float64}`: Internal force density of each point
- `b_ext::Matrix{Float64}`: External force density of each point
- `damage::Vector{Float64}`: Damage of each point
- `n_active_bonds::Vector{Int}`: Number of intact bonds of each point
- `stress::Matrix{Float64}`: Stress tensor of each point
- `von_mises_stress::Vector{Float64}`: Von Mises stress of each point
"""
struct RKCMaterial{CM,K,DM} <: AbstractRKCMaterial{CM,NoCorrection}
    kernel::K
    constitutive_model::CM
    dmgmodel::DM
    maxdmg::Float64
    accuracy_order::Int
    function RKCMaterial(kernel::K, cm::CM, dmgmodel::DM, maxdmg::Real,
                         accuracy_order::Int) where {CM,K,DM}
        return new{CM,K,DM}(kernel, cm, dmgmodel, maxdmg, accuracy_order)
    end
end

function RKCMaterial(; kernel::Function=cubic_b_spline_kernel,
                     model::AbstractConstitutiveModel=LinearElastic(),
                     dmgmodel::AbstractDamageModel=CriticalStretch(), maxdmg::Real=0.85,
                     accuracy_order::Int=1)
    # if kernel !== cubic_b_spline_kernel
    #     msg = "The kernel $kernel is not recommended for use with RKCMaterial."
    #     msg *= "Use the `cubic_b_spline_kernel` function instead.\n"
    #     throw(ArgumentError(msg))
    # end
    if !(accuracy_order in (1, 2))
        msg = "RK kernel only implemented for `accuracy_order ∈ {1,2}`!\n"
        throw(ArgumentError(msg))
    end
    return RKCMaterial(kernel, model, dmgmodel, maxdmg, accuracy_order)
end

function Base.show(io::IO, @nospecialize(mat::RKCMaterial))
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat, (:maxdmg,)))
    return nothing
end

function log_material_property(::Val{:accuracy_order}, mat; indentation)
    return msg_qty("gradient order of accuracy", mat.accuracy_order; indentation)
end

@params RKCMaterial StandardPointParameters

@storage RKCMaterial struct RKCStorage
    @lthfield position::Matrix{Float64}
    @pointfield displacement::Matrix{Float64}
    @pointfield velocity::Matrix{Float64}
    @pointfield velocity_half::Matrix{Float64}
    @pointfield velocity_half_old::Matrix{Float64}
    @pointfield acceleration::Matrix{Float64}
    @htlfield b_int::Matrix{Float64}
    @pointfield b_int_old::Matrix{Float64}
    @pointfield b_ext::Matrix{Float64}
    @pointfield density_matrix::Matrix{Float64}
    @pointfield damage::Vector{Float64}
    bond_active::Vector{Bool}
    @pointfield n_active_bonds::Vector{Int}
    @pointfield update_gradients::Vector{Bool}
    @pointfield stress::Matrix{Float64}
    @pointfield von_mises_stress::Vector{Float64}
    @lthfield defgrad::Matrix{Float64}
    @lthfield weighted_volume::Vector{Float64}
    gradient_weight::Matrix{Float64}
    first_piola_kirchhoff::Matrix{Float64}
end

# function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
#                     ::Val{:displacement})
#     return zeros(3, get_n_points(system))
# end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:update_gradients})
    return ones(Bool, get_n_loc_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:stress})
    return zeros(9, get_n_loc_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:von_mises_stress})
    return zeros(get_n_loc_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:defgrad})
    return zeros(9, get_n_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:weighted_volume})
    return zeros(get_n_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:gradient_weight})
    return zeros(3, get_n_bonds(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:first_piola_kirchhoff})
    return zeros(9, get_n_bonds(system))
end

function initialize!(chunk::BodyChunk{<:BondSystem,<:AbstractRKCMaterial})
    (; system, mat, paramsetup, storage) = chunk
    storage.update_gradients .= true
    for i in each_point_idx(system)
        calc_weights_and_defgrad!(storage, system, mat, paramsetup, 0.0, 0.0, i)
    end
    return nothing
end

function calc_force_density!(dh::ThreadsBodyDataHandler{<:BondSystem,<:AbstractRKCMaterial},
                             t, Δt)
    (; chunks) = dh
    after_fields = rkc_lth_after_fields(chunks[1].mat)
    pre_fields = filter(x -> !in(x, after_fields), loc_to_halo_fields(chunks[1].storage))
    @threads :static for chunk_id in eachindex(chunks)
        exchange_loc_to_halo!(dh, chunk_id, pre_fields)
    end
    @threads :static for chunk_id in eachindex(chunks)
        calc_weights_and_defgrad!(chunks[chunk_id], t, Δt)
    end
    @threads :static for chunk_id in eachindex(chunks)
        exchange_loc_to_halo!(dh, chunk_id, after_fields)
    end
    @threads :static for chunk_id in eachindex(chunks)
        calc_force_density!(chunks[chunk_id], t, Δt)
    end
    @threads :static for chunk_id in eachindex(chunks)
        exchange_halo_to_loc!(dh, chunk_id)
    end
    return nothing
end

rkc_lth_after_fields(::RKCMaterial) = (:defgrad, :weighted_volume)

function calc_force_density!(dh::MPIBodyDataHandler{<:BondSystem,<:AbstractRKCMaterial}, t,
                             Δt)
    (; chunk) = dh
    after_fields = rkc_lth_after_fields(chunk.mat)
    pre_fields = filter(x -> !in(x, after_fields), loc_to_halo_fields(chunk.storage))
    exchange_loc_to_halo!(dh, pre_fields)
    calc_weights_and_defgrad!(chunk, t, Δt)
    exchange_loc_to_halo!(dh, after_fields)
    calc_force_density!(chunk, t, Δt)
    exchange_halo_to_loc!(dh)
    return nothing
end

function calc_weights_and_defgrad!(chunk::BodyChunk{<:BondSystem,<:AbstractRKCMaterial}, t,
                                   Δt)
    (; system, mat, paramsetup, storage) = chunk
    (; dmgmodel) = mat
    storage.n_active_bonds .= 0
    for i in each_point_idx(system)
        calc_failure!(storage, system, mat, dmgmodel, paramsetup, i)
        calc_damage!(storage, system, mat, dmgmodel, paramsetup, i)
        calc_weights_and_defgrad!(storage, system, mat, paramsetup, t, Δt, i)
    end
    return nothing
end

function calc_damage!(storage::AbstractStorage, system::AbstractBondSystem,
                      mat::AbstractRKCMaterial, dmgmodel::AbstractDamageModel,
                      paramsetup::AbstractParameterSetup, i)
    (; n_neighbors) = system
    (; n_active_bonds, damage, update_gradients) = storage
    old_damage = damage[i]
    new_damage = 1 - n_active_bonds[i] / n_neighbors[i]
    if new_damage > old_damage
        update_gradients[i] = true
    else
        update_gradients[i] = false
    end
    damage[i] = new_damage
    return nothing
end

function update_gradients(::AbstractRKCMaterial, storage::AbstractStorage, i)
    return storage.update_gradients[i]
end

function calc_weights_and_defgrad!(storage::AbstractStorage, system::AbstractBondSystem,
                                   mat::AbstractRKCMaterial,
                                   paramsetup::AbstractParameterSetup, t, Δt, i)
    params = get_params(paramsetup, i)
    if update_gradients(mat, storage, i)
        rkc_weights!(storage, system, mat, params, t, Δt, i)
    end
    rkc_defgrad!(storage, system, mat, params, t, Δt, i)
    return nothing
end

function rkc_weights!(storage::AbstractStorage, system::AbstractBondSystem,
                     mat::AbstractRKCMaterial, params::AbstractPointParameters, t, Δt, i)
    (; bonds, volume) = system
    (; bond_active, gradient_weight, weighted_volume, update_gradients) = storage
    (; accuracy_order) = mat
    q_dim = get_q_dim(accuracy_order)

    Q∇ᵀ = get_q_triangle(accuracy_order)

    (; δ) = params

    # calculate moment matrix M
    M = zero(SMatrix{q_dim,q_dim,Float64,q_dim*q_dim})
    wi = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        # Q = get_monomial_vector(accuracy_order, ΔXij, δ)
        Q = get_monomial_vector(accuracy_order, ΔXij, 1.0)
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        temp = ωij * volume[j]
        M += temp * (Q * Q')
        wi += temp
    end
    weighted_volume[i] = wi

    # calculate inverse of moment matrix, must be a full rank matrix!
    threshold = 1e-9 * δ * δ * δ
    Minv = invreg(M, threshold)
    # Minv = inv(M)

    # calculate gradient weights Φ
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        # Q = get_monomial_vector(accuracy_order, ΔXij, δ)
        Q = get_monomial_vector(accuracy_order, ΔXij, 1.0)
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        # temp = ωij / δ
        # temp = ωij / δ * volume[j]
        temp = ωij * volume[j]
        MinvQ = Minv * Q
        Φ = temp * (Q∇ᵀ * MinvQ)
        update_vector!(gradient_weight, bond_id, Φ)
    end

    # gradients are evaluated and do not need to be updated anymore
    update_gradients[i] = false

    return nothing
end

# function my_rkc_weights!(storage::AbstractStorage, system::AbstractBondSystem,
#                          mat::AbstractRKCMaterial, params::AbstractPointParameters, t, Δt, i)
#     (; bonds, volume) = system
#     (; bond_active, gradient_weight, weighted_volume, update_gradients) = storage
#     (; accuracy_order) = mat

#     q_dim = get_q_dim(accuracy_order)
#     Q∇ᵀ = get_q_triangle(accuracy_order)

#     (; δ) = params

#     # calculate moment matrix M
#     # M2 = zero(MMatrix{q_dim,q_dim,Float64,q_dim*q_dim})
#     M = zero(SMatrix{q_dim,q_dim,Float64,q_dim*q_dim})
#     wi = 0.0
#     for bond_id in each_bond_idx(system, i)
#         bond = bonds[bond_id]
#         j = bond.neighbor
#         ΔXij = get_vector_diff(system.position, i, j)
#         Q = get_monomial_vector(accuracy_order, ΔXij, δ)
#         ωij = kernel(system, bond_id) * bond_active[bond_id]
#         temp = ωij * volume[j]
#         # for jj in 1:q_dim
#         #     for ii in 1:q_dim
#         #         M2[ii, jj] += temp * Q[ii] * Q[jj]
#         #     end
#         # end
#         M += temp * (Q * Q')
#         wi += temp
#     end
#     weighted_volume[i] = wi

#     # @show M ≈ M2

#     # error("stop")

#     # calculate inverse of moment matrix, must be a full rank matrix!
#     threshold = 1e-6 * δ * δ * δ
#     Minv = invreg(M, threshold)

#     # calculate gradient weights Φ
#     for bond_id in each_bond_idx(system, i)
#         bond = bonds[bond_id]
#         j = bond.neighbor
#         ΔXij = get_vector_diff(system.position, i, j)
#         Q = get_monomial_vector(accuracy_order, ΔXij, δ)
#         ωij = kernel(system, bond_id) * bond_active[bond_id]
#         temp = ωij / δ
#         # temp = ωij / δ * volume[j]
#         MinvQ = Minv * Q
#         Φ = temp * (Q∇ᵀ * MinvQ)
#         Φ1 = 0.0
#         Φ2 = 0.0
#         Φ3 = 0.0
#         temp = ωij
#         for jj in 1:q_dim
#             Φ1 += temp * Minv[jj, 1] * Q[jj] / δ
#             Φ2 += temp * Minv[jj, 2] * Q[jj] / δ
#             Φ3 += temp * Minv[jj, 3] * Q[jj] / δ
#         end
#         Φ_2 = SVector{3,Float64}(Φ1, Φ2, Φ3)
#         @show Φ ≈ Φ_2
#         update_vector!(gradient_weight, bond_id, Φ)
#     end

#     # gradients are evaluated and do not need to be updated anymore
#     update_gradients[i] = false

#     return nothing
# end

# @inline each_monomial(N::Int) = each_monomial(Val(N))

# @inline function each_monomial(::Val{N}) where {N}
#     monomials = Tuple((p1,p2,p3) for i in 1:N
#                                  for p1 in i:-1:0
#                                  for p2 in (i - p1):-1:0
#                                  for p3 in i - p1 - p2)
#     return ((0, 0, 0), monomials...)
# end

# @inline function each_monomial(::Val{1})
#     return ((1, 0, 0), (0, 1, 0), (0, 0, 1))
# end

# @inline function each_monomial(::Val{2})
#     # return ((1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0),
#     #         (0, 1, 1), (0, 0, 2))
#     return ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), (0, 2, 0), (0, 0, 2))
# end

# @inline get_q_dim(accuracy_order) = length(each_monomial(accuracy_order))
@inline get_q_dim(N::Int) = get_q_dim(Val(N))
@inline get_q_dim(::Val{1}) = 3
@inline get_q_dim(::Val{2}) = 9

@inline function get_monomial_vector(N::Int, ΔX::AbstractArray, δ::Real)
    return get_monomial_vector(Val(N), ΔX, δ)
end

# @inline function get_monomial_vector(n::Val{N}, ΔX::AbstractArray, δ::Real) where {N}
#     a1, a2, a3 = ΔX[1] / δ, ΔX[2] / δ, ΔX[3] / δ
#     _Q = (a1^p1 * a2^p2 * a3^p3 for (p1, p2, p3) in each_monomial(n))
#     q_dims = get_q_dim(n)
#     Q = SVector{q_dims,eltype(ΔX)}(_Q...)
#     return Q
# end

# @inline function get_monomial_vector(n::Val{N}, ΔX::AbstractArray, δ::Real) where {N}
#     a1, a2, a3 = ΔX[1] / δ, ΔX[2] / δ, ΔX[3] / δ
#     _Q = (a1^p1 * a2^p2 * a3^p3 for (p1, p2, p3) in each_monomial(n))
#     q_dims = get_q_dim(n)
#     Q = SVector{q_dims,eltype(ΔX)}(_Q...)
#     return Q
# end

@inline function get_monomial_vector(::Val{1}, ΔX::AbstractArray, δ::Real)
    a1, a2, a3 = ΔX[1] / δ, ΔX[2] / δ, ΔX[3] / δ
    Q = SVector{3,eltype(ΔX)}(a1, a2, a3)
    return Q
end

@inline function get_monomial_vector(::Val{2}, ΔX::AbstractArray, δ::Real)
    a1, a2, a3 = ΔX[1] / δ, ΔX[2] / δ, ΔX[3] / δ
    Q = SVector{9,eltype(ΔX)}(a1, a2, a3, a1*a1, a1*a2, a1*a3, a2*a2, a2*a3, a3*a3)
    # Q = SVector{10,eltype(ΔX)}(1.0, a1, a2, a3, a1*a1, a1*a2, a1*a3, a2*a2, a2*a3, a3*a3)
    # Q = SVector{7,eltype(ΔX)}(1.0, a1, a2, a3, a1*a1, a2*a2, a3*a3)
    return Q
end

@inline get_q_triangle(N::Int) = get_q_triangle(Val(N))

function get_q_triangle(::Val{1})
    Q∇ᵀ = SMatrix{3,3,Int,9}(1, 0, 0,
                             0, 1, 0,
                             0, 0, 1)
    return Q∇ᵀ
end

function get_q_triangle(::Val{2})
    Q∇ᵀ = SMatrix{3,9,Int,27}(1, 0, 0,
                              0, 1, 0,
                              0, 0, 1,
                              0, 0, 0,
                              0, 0, 0,
                              0, 0, 0,
                              0, 0, 0,
                              0, 0, 0,
                              0, 0, 0)
    return Q∇ᵀ
end

function get_q_triangle(::Val{N}) where {N}
    msg = "RK kernel only implemented for `accuracy_order ∈ {1,2}`!\n"
    return throw(ArgumentError(msg))
end

function rkc_defgrad!(storage::AbstractStorage, system::AbstractBondSystem,
                      mat::AbstractRKCMaterial, params::AbstractPointParameters, t, Δt, i)
    (; bonds, volume) = system
    (; bond_active, defgrad, gradient_weight) = storage

    F = SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    # F2 = MMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        Δxij = get_vector_diff(storage.position, i, j)
        # ΔUij = get_vector_diff(storage.displacement, i, j)
        ΔUij = Δxij - ΔXij
        Φij = get_vector(gradient_weight, bond_id)
        F += (ΔUij * Φij')# * volume[j]
        # for ii in 1:3
        #     F2[ii, 1] += ΔUij[ii] * Φij[1] * volume[j]
        #     F2[ii, 2] += ΔUij[ii] * Φij[2] * volume[j]
        #     F2[ii, 3] += ΔUij[ii] * Φij[3] * volume[j]
        # end
    end

    # (; accuracy_order) = mat
    # _F = SMatrix{3,3,Float64,9}(0.0, 0.0, 0.0,
    #                             0.0, 0.0, 0.0,
    #                             0.0, 0.0, 0.0)
    # q_dim = get_q_dim(accuracy_order)
    # M = zero(SMatrix{q_dim,q_dim,Float64,q_dim*q_dim})
    # Q∇ᵀ = get_q_triangle(accuracy_order)
    # for bond_id in each_bond_idx(system, i)
    #     bond = bonds[bond_id]
    #     j = bond.neighbor
    #     ΔXij = get_vector_diff(system.position, i, j)
    #     Δxij = get_vector_diff(storage.position, i, j)
    #     # ΔUij = get_vector_diff(storage.displacement, i, j)
    #     # ΔUij = Δxij - ΔXij
    #     # Φij = get_vector(gradient_weight, bond_id)
    #     # _F += (ΔUij * Φij') * volume[j]

    #     Q = get_monomial_vector(accuracy_order, ΔXij, 1.0)
    #     @assert Q ≈ ΔXij
    #     ωij = kernel(system, bond_id) * bond_active[bond_id]
    #     temp = ωij * volume[j]
    #     M += temp * (Q * Q')
    #     _F += temp * (Δxij * Q')
    # end
    # δ = params.δ
    # threshold = 1e-6 * δ * δ * δ
    # Minv = invreg(M, threshold)
    # Minv = inv(M)
    # F = _F * Minv

    # (; bonds, volume) = system
    # (; bond_active, defgrad, weighted_volume) = storage

    # K = zero(SMatrix{3,3,Float64,9})
    # _F = zero(SMatrix{3,3,Float64,9})
    # wi = 0.0
    # for bond_id in each_bond_idx(system, i)
    #     bond = bonds[bond_id]
    #     j = bond.neighbor
    #     ΔXij = get_vector_diff(system.position, i, j)
    #     Δxij = get_vector_diff(storage.position, i, j)
    #     ωij = kernel(system, bond_id) * bond_active[bond_id]
    #     temp = ωij * volume[j]
    #     K += temp * (ΔXij * ΔXij')
    #     _F += temp * (Δxij * ΔXij')
    #     wi += temp
    # end
    # Kinv = inv(K)
    # F = _F * Kinv


    update_tensor!(defgrad, i, F)

    return nothing
end

function calc_force_density!(chunk::BodyChunk{<:BondSystem,<:AbstractRKCMaterial}, t, Δt)
    (; system, mat, paramsetup, storage) = chunk
    storage.b_int .= 0
    for i in each_point_idx(chunk)
        force_density_point!(storage, system, mat, paramsetup, t, Δt, i)
    end
    nancheck(chunk, t)
    return nothing
end

function force_density_point!(storage::AbstractStorage, system::BondSystem,
                              mat::AbstractRKCMaterial, params::AbstractPointParameters, t,
                              Δt, i)
    ∑P = rkc_stress_integral!(storage, system, mat, params, t, Δt, i)
    rkc_force_density!(storage, system, mat, params, ∑P, t, Δt, i)
    return nothing
end

# function force_density_point!(storage::AbstractStorage, system::BondSystem,
#                               mat::AbstractRKCMaterial, params::AbstractPointParameters, t,
#                               Δt, i)
#     (; bonds, volume) = system
#     (; bond_active, gradient_weight, b_int) = storage
#     F = get_tensor(storage.defgrad, i)
#     P = first_piola_kirchhoff(mat.constitutive_model, storage, params, F)
#     for bond_id in each_bond_idx(system, i)
#         bond = bonds[bond_id]
#         j = bond.neighbor
#         if bond_active[bond_id]
#             Ψ = get_vector(gradient_weight, bond_id)
#             tij = P * Ψ
#             update_add_vector!(b_int, i, tij * volume[j])
#             update_add_vector!(b_int, j, -tij * volume[i])
#         end
#     end
#     return nothing
# end

function rkc_stress_integral!(storage::AbstractStorage, system::AbstractBondSystem,
                              mat::AbstractRKCMaterial, params::AbstractPointParameters, t,
                              Δt, i)
    (; bonds) = system
    (; bond_active, defgrad, weighted_volume) = storage
    Fi = get_tensor(defgrad, i)
    # too_much_damage!(storage, system, mat, Fi, i) && return nothing
    wi = weighted_volume[i]
    ∑P = SMatrix{3,3,Float64,9}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length
            ΔXij = get_vector_diff(system.position, i, j)
            Δxij = get_vector_diff(storage.position, i, j)
            Fj = get_tensor(defgrad, j)
            Fb = 0.5 * (Fi + Fj)
            ΔXijLL = ΔXij' / (L * L)
            ΔFij = (Δxij - Fb * ΔXij) * ΔXijLL
            Fij = Fb + ΔFij
            Pij = calc_first_piola_kirchhoff!(storage, mat, params, Fij, bond_id)
            Tempij = I - ΔXij * ΔXijLL
            wj = weighted_volume[j]
            ϕ = (wi > 0 && wj > 0) ? (0.5 / wi + 0.5 / wj) : 0.0
            ω̃ij = kernel(system, bond_id) * ϕ
            ∑Pij = ω̃ij * (Pij * Tempij)
            ∑P += ∑Pij
        end
    end
    return ∑P
end

function rkc_force_density!(storage::AbstractStorage, system::AbstractBondSystem,
                            mat::AbstractRKCMaterial, params::AbstractPointParameters,
                            ∑P, t, Δt, i)
    (; bonds, volume) = system
    (; bond_active, gradient_weight, weighted_volume) = storage
    wi = weighted_volume[i]
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length
            ΔXij = get_vector_diff(system.position, i, j)
            Pij = get_tensor(storage.first_piola_kirchhoff, bond_id)
            Φij = get_vector(gradient_weight, bond_id)
            wj = weighted_volume[j]
            ϕ = (wi > 0 && wj > 0) ? (0.5 / wi + 0.5 / wj) : 0.0
            ω̃ij = kernel(system, bond_id) * ϕ
            tij = ω̃ij / (L * L) * (Pij * ΔXij) + ∑P * Φij
            update_add_vector!(storage.b_int, i, tij * volume[j])
            update_add_vector!(storage.b_int, j, -tij * volume[i])
        end
    end
    return nothing
end

function calc_first_piola_kirchhoff!(storage::RKCStorage, mat::RKCMaterial,
                                     params::StandardPointParameters, F, bond_id)
    P = first_piola_kirchhoff(mat.constitutive_model, storage, params, F)
    update_tensor!(storage.first_piola_kirchhoff, bond_id, P)
    return P
end

# function too_much_damage!(storage::AbstractStorage, system::BondSystem,
#                           mat::AbstractRKCMaterial, F, i)
#     if storage.damage[i] > mat.maxdmg || containsnan(F)
#         # kill all bonds of this point
#         storage.bond_active[each_bond_idx(system, i)] .= false
#         return true
#     end
#     return false
# end
