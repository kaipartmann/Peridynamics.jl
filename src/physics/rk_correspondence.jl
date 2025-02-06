"""
    RKCMaterial(; maxdmg, zem)

A material type used to assign the material of a [`Body`](@ref) with the local continuum
consistent (correspondence) formulation of non-ordinary state-based peridynamics.

# Keywords
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. If this value is
    exceeded, all bonds of that point are broken because the deformation gradient would then
    possibly contain `NaN` values.
    (default: `0.95`)
- `zem::AbstractZEMStabilization`: zero-energy mode stabilization. The
    stabilization algorithm of Silling (2017) is used as default.
    (default: `ZEMSilling()`)

!!! note "Stability of fracture simulations"
    This formulation is known to be not suitable for fracture simultations without
    stabilization of the zero-energy modes. Therefore be careful when doing fracture
    simulations and try out different parameters for `maxdmg` and `zem`.

# Examples

```julia-repl
julia> mat = RKCMaterial()
RKCMaterial(maxdmg=0.95, zem_fac=ZEMSilling())
```

---

```julia
RKCMaterial
```

Material type for the local continuum consistent (correspondence) formulation of
non-ordinary state-based peridynamics.

# Fields
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. See the
    constructor docs for more informations.
- `zem_fac::Float64`: Correction factor used for zero-energy mode stabilization. See the
    constructor docs for more informations.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `RKCMaterial`, then the following
parameters are allowed:
- `horizon::Float64`: Radius of point interactions
- `rho::Float64`: Density
- `E::Float64`: Young's modulus
- `nu::Float64`: Poisson's ratio
- `Gc::Float64`: Critical energy release rate
- `epsilon_c::Float64`: Critical strain

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
"""
struct RKCMaterial{CM,K} <: AbstractCorrespondenceMaterial{CM,NoCorrection}
    kernel::K
    constitutive_model::CM
    maxdmg::Float64
    accuracy_order::Int
    function RKCMaterial(kernel::K, cm::CM, maxdmg::Real, accuracy_order::Int) where {CM,K}
        return new{CM,K}(kernel, cm, maxdmg, accuracy_order)
    end
end

function Base.show(io::IO, @nospecialize(mat::RKCMaterial))
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat, (:maxdmg,)))
    return nothing
end

function RKCMaterial(; kernel::Function=cubic_b_spline,
                     model::AbstractConstitutiveModel=SaintVenantKirchhoff(),
                     maxdmg::Real=0.85, accuracy_order::Int=2)
    return RKCMaterial(kernel, model, maxdmg, accuracy_order)
end

struct RKCPointParameters <: AbstractPointParameters
    δ::Float64
    rho::Float64
    E::Float64
    nu::Float64
    G::Float64
    K::Float64
    λ::Float64
    μ::Float64
    Gc::Float64
    εc::Float64
    bc::Float64
end

function RKCPointParameters(mat::RKCMaterial, p::Dict{Symbol,Any})
    (; δ, rho, E, nu, G, K, λ, μ) = get_required_point_parameters(mat, p)
    (; Gc, εc) = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return RKCPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params RKCMaterial RKCPointParameters

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
    @pointfield damage_changed::Vector{Bool}
    @pointfield pk1::Matrix{Float64}
    # @pointfield stress::Matrix{Float64}
    # @pointfield von_mises_stress::Vector{Float64}
    @lthfield defgrad::Matrix{Float64}
    @lthfield weighted_volume::Vector{Float64}
    gradient_weight::Matrix{Float64}
    first_piola_kirchhoff::Matrix{Float64}
end

# function init_field(::RKCMaterial, ::AbstractTimeSolver, system::BondSystem,
#                     ::Val{:displacement})
#     return zeros(3, get_n_points(system))
# end

function init_field(::RKCMaterial, ::AbstractTimeSolver, system::BondSystem, ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function init_field(::RKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:damage_changed})
    return ones(Bool, get_n_loc_points(system))
end

function init_field(::RKCMaterial, ::AbstractTimeSolver, system::BondSystem, ::Val{:pk1})
    return zeros(9, get_n_loc_points(system))
end

# function init_field(::RKCMaterial, ::AbstractTimeSolver, system::BondSystem,
#                     ::Val{:von_mises_stress})
#     return zeros(get_n_loc_points(system))
# end

function init_field(::RKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:defgrad})
    return zeros(9, get_n_points(system))
end

function init_field(::RKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:weighted_volume})
    return zeros(get_n_points(system))
end

function init_field(::RKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:gradient_weight})
    return zeros(3, get_n_bonds(system))
end

function init_field(::RKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:first_piola_kirchhoff})
    return zeros(9, get_n_bonds(system))
end

function initialize!(chunk::BodyChunk{<:BondSystem,<:RKCMaterial})
    chunk.storage.damage_changed .= true
    return nothing
end

# @inline function calc_damage!(chunk::BodyChunk{<:BondSystem,<:RKCMaterial})
#     (; n_neighbors) = chunk.system
#     (; n_active_bonds, damage, damage_changed) = chunk.storage
#     for point_id in each_point_idx(chunk)
#         old_damage = damage[point_id]
#         new_damage = 1 - n_active_bonds[point_id] / n_neighbors[point_id]
#         if new_damage > old_damage
#             damage_changed[point_id] = true
#         else
#             damage_changed[point_id] = false
#         end
#         damage[point_id] = new_damage
#     end
#     return nothing
# end

@inline function calc_damage!(chunk::BodyChunk{<:BondSystem,<:RKCMaterial})
    (; system, storage, paramsetup) = chunk
    (; n_neighbors, bonds) = system
    (; n_active_bonds, damage, damage_changed) = storage
    for i in each_point_idx(chunk)
        # params = get_params(paramsetup, i)
        # for bond_id in each_bond_idx(system, i)
        #     bond = bonds[bond_id]
        #     j, L = bond.neighbor, bond.length
        #     Δxij = get_coordinates_diff(storage, i, j)
        #     l = norm(Δxij)
        #     ε = (l - L) / L
        #     stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)
        # end
        old_damage = damage[i]
        new_damage = 1 - n_active_bonds[i] / n_neighbors[i]
        if new_damage > old_damage
            damage_changed[i] = true
        else
            damage_changed[i] = false
        end
        damage[i] = new_damage
    end
    return nothing
end

function calc_force_density!(dh::ThreadsBodyDataHandler{<:BondSystem,<:RKCMaterial}, t, Δt)
    @threads :static for chunk_id in eachindex(dh.chunks)
        exchange_loc_to_halo!(dh, chunk_id, :position)
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        calc_weights_and_defgrad!(dh.chunks[chunk_id], t, Δt)
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        exchange_loc_to_halo!(dh, chunk_id, (:defgrad, :weighted_volume))
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        calc_force_density!(chunk, t, Δt)
        nancheck(chunk, t)
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        exchange_halo_to_loc!(dh, chunk_id)
    end
    return nothing
end

function calc_force_density!(dh::MPIBodyDataHandler{<:BondSystem,<:RKCMaterial}, t, Δt)
    (; chunk) = dh
    exchange_loc_to_halo!(dh, (:position,))
    calc_weights_and_defgrad!(chunk, t, Δt)
    exchange_loc_to_halo!(dh, (:defgrad, :weighted_volume))
    calc_force_density!(chunk, t, Δt)
    nancheck(chunk, t)
    exchange_halo_to_loc!(dh)
    return nothing
end

function calc_weights_and_defgrad!(chunk::BodyChunk{<:BondSystem,<:RKCMaterial}, t, Δt)
    (; system, mat, paramsetup, storage) = chunk
    calc_gradient_weights!(storage, system, mat, paramsetup)
    calc_deformation_gradients!(storage, system, mat, paramsetup, t, Δt)
    return nothing
end

function calc_gradient_weights!(storage::RKCStorage, system::BondSystem, mat::RKCMaterial,
                                paramsetup::AbstractParameterSetup)
    (; bonds, volume) = system
    (; bond_active, gradient_weight, damage_changed, weighted_volume) = storage
    (; accuracy_order) = mat
    q_dim = get_q_dim(accuracy_order)

    Q∇ᵀ = get_q_triangle(accuracy_order)

    for i in each_point_idx(system)
        if damage_changed[i]
            params = get_params(paramsetup, i)
            (; δ) = params

            # calculate moment matrix M
            M = zero(SMatrix{q_dim,q_dim,Float64,q_dim*q_dim})
            wi = 0.0
            for bond_id in each_bond_idx(system, i)
                bond = bonds[bond_id]
                j = bond.neighbor
                ΔXij = get_diff(system.position, i, j)
                Q = get_monomial_vector(accuracy_order, ΔXij, δ)
                ωij = kernel(system, bond_id) #* bond_active[bond_id]
                temp = ωij * volume[j]
                M += temp * (Q * Q')
                wi += temp
            end
            weighted_volume[i] = wi

            # calculate inverse of moment matrix, must be a full rank matrix!
            threshold = 1e-6 * δ * δ * δ
            Minv = invreg(M, threshold)

            # calculate gradient weights Φ
            for bond_id in each_bond_idx(system, i)
                bond = bonds[bond_id]
                j = bond.neighbor
                ΔXij = get_diff(system.position, i, j)
                Q = get_monomial_vector(accuracy_order, ΔXij, δ)
                ωij = kernel(system, bond_id) * bond_active[bond_id]
                temp = ωij / δ * volume[j]
                MinvQ = Minv * Q
                Φ = temp * (Q∇ᵀ * MinvQ)
                #---
                # Φ1 = 0.0
                # Φ2 = 0.0
                # Φ3 = 0.0
                # for m in 1:q_dim
                #     Φ1 += temp * Minv[1, m] * Q[m]
                #     Φ2 += temp * Minv[2, m] * Q[m]
                #     Φ3 += temp * Minv[3, m] * Q[m]
                # end
                # Φm = SVector{3,Float64}(Φ1, Φ2, Φ3)
                # @assert isapprox(Φ, Φm)
                # @show Φ, Φm
                #---
                update_vector!(gradient_weight, bond_id, Φ)
            end
        end
    end

    return nothing
end

@inline each_monomial(N::Int) = each_monomial(Val(N))

@inline function each_monomial(::Val{N}) where {N}
    monomials = Tuple((p1,p2,p3) for i in 1:N
                                 for p1 in i:-1:0
                                 for p2 in (i - p1):-1:0
                                 for p3 in i - p1 - p2)
    return monomials
end

@inline function each_monomial(::Val{1})
    return ((1, 0, 0), (0, 1, 0), (0, 0, 1))
end

@inline function each_monomial(::Val{2})
    return ((1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0),
            (0, 1, 1), (0, 0, 2))
end

@inline get_q_dim(accuracy_order) = length(each_monomial(accuracy_order))

@inline function get_monomial_vector(N::Int, ΔX::AbstractArray, δ::Real)
    return get_monomial_vector(Val(N), ΔX, δ)
end

@inline function get_monomial_vector(n::Val{N}, ΔX::AbstractArray, δ::Real) where {N}
    a1, a2, a3 = ΔX[1] / δ, ΔX[2] / δ, ΔX[3] / δ
    _Q = (a1^p1 * a2^p2 * a3^p3 for (p1, p2, p3) in each_monomial(n))
    q_dims = get_q_dim(n)
    Q = SVector{q_dims,eltype(ΔX)}(_Q...)
    return Q
end

@inline function get_monomial_vector(::Val{1}, ΔX::AbstractArray, δ::Real)
    a1, a2, a3 = ΔX[1] / δ, ΔX[2] / δ, ΔX[3] / δ
    Q = SVector{3,eltype(ΔX)}(a1, a2, a3)
    return Q
end

@inline function get_monomial_vector(::Val{2}, ΔX::AbstractArray, δ::Real)
    a1, a2, a3 = ΔX[1] / δ, ΔX[2] / δ, ΔX[3] / δ
    Q = SVector{9,eltype(ΔX)}(a1, a2, a3, a1*a1, a1*a2, a1*a3, a2*a2, a2*a3, a3*a3)
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

# function invert_moment_matrix(M::StaticMatrix{N,N,T}, threshold::Real) where {N,T}
#     U, S, V = svd(M)
#     Sinvreg = SVector{N,T}((s > threshold ? 1/s : 0) for s in S)
#     Sinv = Diagonal{T,SVector{N,T}}(Sinvreg)
#     Minv = V * Sinv * U'
#     return Minv
# end

function calc_deformation_gradients!(storage::RKCStorage, system::BondSystem, ::RKCMaterial,
                                     ::AbstractParameterSetup, t, Δt)
    (; bonds, volume) = system
    (; gradient_weight, defgrad, bond_active) = storage

    for i in each_point_idx(system)
        F = SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
        for bond_id in each_bond_idx(system, i)
            # if bond_active[bond_id]
                bond = bonds[bond_id]
                j = bond.neighbor
                ΔXij = get_vector_diff(system.position, i, j)
                Δxij = get_vector_diff(storage.position, i, j)
                ΔUij = Δxij - ΔXij
                Φij = get_vector(gradient_weight, bond_id)
                F += (ΔUij * Φij') * volume[j]
            # end
        end
        update_tensor!(defgrad, i, F)
    end

    return nothing
end

function force_density_point!(storage::RKCStorage, system::BondSystem, mat::RKCMaterial,
                              paramhandler::AbstractParameterHandler, t, Δt, i)
    params = get_params(paramhandler, i)
    force_density_point!(storage, system, mat, params, t, Δt, i)
    return nothing
end

function force_density_point!(storage::RKCStorage, system::BondSystem, mat::RKCMaterial,
                              params::RKCPointParameters, t, Δt, i)
    # (; bonds, volume) = system
    # (; gradient_weight, defgrad, bond_active, weighted_volume) = storage
    # (; first_piola_kirchhoff) = storage
    # Fi = get_tensor(defgrad, i)
    # too_much_damage!(storage, system, mat, Fi, i) && return nothing
    # Pi = calc_first_piola_kirchhoff(storage, mat, params, Fi)

    # Standard RK
    # for bond_id in each_bond_idx(system, i)
    #     bond = bonds[bond_id]
    #     j, L = bond.neighbor, bond.length
    #     Δxij = get_coordinates_diff(storage, i, j)
    #     l = norm(Δxij)
    #     ε = (l - L) / L
    #     stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)

    #     if bond_active[bond_id]
    #         Fj = get_tensor(defgrad, j)
    #         Fb = 0.5 * (Fi + Fj)
    #         ΔXij = get_coordinates_diff(system, i, j)
    #         ΔFij = (Δxij - Fb * ΔXij) * ΔXij' / (L * L)
    #         Fij = Fb + ΔFij
    #         Pij = calc_first_piola_kirchhoff(storage, mat, params, Fij)
    #         Φij = get_vector(gradient_weight, bond_id)
    #         tij = Pij * Φij
    #         update_add_b_int!(storage, i, tij * volume[j])
    #         update_add_b_int!(storage, j, -tij * volume[i])
    #     end
    # end

    # BA-Stabilized RK
    # for bond_id in each_bond_idx(system, i)
    #     bond = bonds[bond_id]
    #     j, L = bond.neighbor, bond.length
    #     ΔXij = get_coordinates_diff(system, i, j)
    #     Δxij = get_coordinates_diff(storage, i, j)
    #     l = norm(Δxij)
    #     ε = (l - L) / L
    #     stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)

    #     Φij = get_vector(gradient_weight, bond_id)
    #     Fj = get_tensor(defgrad, j)
    #     Fb = 0.5 * (Fi + Fj)
    #     ΔFij = (Δxij - Fb * ΔXij) * ΔXij' / (L * L)
    #     Fij = Fj + ΔFij
    #     Pij = calc_first_piola_kirchhoff(storage, mat, params, Fij) * bond_active[bond_id]
    #     tij = (Pij * Φij)
    #     update_add_b_int!(storage, i, tij * volume[j])
    #     update_add_b_int!(storage, j, -(Pi * Φij) * volume[i])
    # end

    # Nodal quadrature peridigm
    # (; bonds, volume) = system
    # (; bond_active, gradient_weight, defgrad, weighted_volume) = storage

    # Fi = get_tensor(defgrad, i)
    # # too_much_damage!(storage, system, mat, Fi, i) && return nothing

    # # Stress integral
    # wi = weighted_volume[i]
    # ∑P = SMatrix{3,3,Float64,9}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # for bond_id in each_bond_idx(system, i)
    #     bond = bonds[bond_id]
    #     j, L = bond.neighbor, bond.length
    #     ΔXij = get_coordinates_diff(system, i, j)
    #     Δxij = get_coordinates_diff(storage, i, j)
    #     l = norm(Δxij)
    #     ε = (l - L) / L
    #     stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)

    #     if bond_active[bond_id]
    #         Fj = get_tensor(defgrad, j)
    #         Fb = 0.5 * (Fi + Fj)
    #         ΔXijLL = ΔXij' / (L * L)
    #         ΔFij = (Δxij - Fb * ΔXij) * ΔXijLL
    #         Fij = Fb + ΔFij
    #         Pij = calc_first_piola_kirchhoff(storage, mat, params, Fij)
    #         update_tensor!(storage.first_piola_kirchhoff, bond_id, Pij)
    #         Tempij = I - ΔXij * ΔXijLL
    #         # ω̃ij = kernel(system, bond_id) * (0.5 / wi + 0.5 / weighted_volume[j])
    #         wj = weighted_volume[j]
    #         ϕ = (wi > 0 && wj > 0) ? (0.5 / wi + 0.5 / wj) : 0.0
    #         ω̃ij = kernel(system, bond_id) * ϕ
    #         ∑Pij = ω̃ij * (Pij * Tempij)
    #         ∑P += ∑Pij
    #     end
    # end

    # for bond_id in each_bond_idx(system, i)
    #     if bond_active[bond_id]
    #         bond = bonds[bond_id]
    #         j, L = bond.neighbor, bond.length
    #         ΔXij = get_coordinates_diff(system, i, j)
    #         Pij = get_tensor(storage.first_piola_kirchhoff, bond_id)
    #         Φij = get_vector(gradient_weight, bond_id)
    #         # ω̃ij = kernel(system, bond_id) * (0.5 / wi + 0.5 / weighted_volume[j])
    #         wj = weighted_volume[j]
    #         ϕ = (wi > 0 && wj > 0) ? (0.5 / wi + 0.5 / wj) : 0.0
    #         ω̃ij = kernel(system, bond_id) * ϕ
    #         tij = ω̃ij / (L * L) * (Pij * ΔXij) + ∑P * Φij
    #         update_add_b_int!(storage, i, tij * volume[j])
    #         update_add_b_int!(storage, j, -tij * volume[i])
    #     end
    # end
    # return nothing


    ## My understanding of the equations
    # wi = weighted_volume[i]
    # Vi = volume[i]
    # ∑ωPV = SMatrix{3,3,Float64,9}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # ∑ωPcorr = SMatrix{3,3,Float64,9}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # for bond_id in each_bond_idx(system, i)
    #     bond = bonds[bond_id]
    #     j, L = bond.neighbor, bond.length
    #     Δxij = get_coordinates_diff(storage, i, j)
    #     l = norm(Δxij)
    #     ε = (l - L) / L
    #     stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)

    #     if bond_active[bond_id]
    #         Vj = volume[j]
    #         ΔXij = get_coordinates_diff(system, i, j)
    #         Fj = get_tensor(defgrad, j)
    #         Fb = 0.5 * (Fi + Fj)
    #         ΔFij = (Δxij - Fb * ΔXij) * ΔXij' / (L * L)
    #         Fij = Fb + ΔFij
    #         Pij = calc_first_piola_kirchhoff(storage, mat, params, Fij)
    #         update_tensor!(first_piola_kirchhoff, bond_id, Pij)

    #         ωij = kernel(system, bond_id)
    #         ωPV = ωij * Pij * Vj
    #         ∑ωPV += ωPV
    #         ∑ωPcorr += ωPV * (ΔXij * ΔXij') / (L * L)
    #         # wj = weighted_volume[j]
    #         # ϕ = volume[j] * ωij * (0.5 / wi + 0.5 / wj)
    #         # SI += ϕ * (Pij * Tempij)
    #     end
    # end

    ### my newest approach
    (; bonds, volume) = system
    (; gradient_weight, defgrad, bond_active, weighted_volume) = storage
    Fi = get_tensor(defgrad, i)
    Pi = calc_first_piola_kirchhoff(storage, mat, params, Fi)

    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        # if bond_active[bond_id]
            Fj = get_tensor(defgrad, j)
            Pj = calc_first_piola_kirchhoff(storage, mat, params, Fj)
            Φij = get_vector(gradient_weight, bond_id)
            tij = (Pj - Pi) * Φij
            # tij = Pi * Φij
            update_add_b_int!(storage, i, tij * volume[j])
            update_add_b_int!(storage, j, -tij * volume[i])
        # end
    end

    return nothing
end

function calc_first_piola_kirchhoff(storage::RKCStorage, mat::RKCMaterial,
                                    params::RKCPointParameters, F)
    P = first_piola_kirchhoff(mat.constitutive_model, storage, params, F)
    return P
end

# function rkc_force_density!(storage::RKCStorage, system::BondSystem, mat::RKCMaterial,
#                             params::RKCPointParameters, PKinv::SMatrix, i)
#     (; bonds, volume) = system
#     (; bond_active) = storage
#     for bond_id in each_bond_idx(system, i)
#         bond = bonds[bond_id]
#         j, L = bond.neighbor, bond.length
#         ΔXij = get_coordinates_diff(system, i, j)
#         Δxij = get_coordinates_diff(storage, i, j)
#         l = norm(Δxij)
#         ε = (l - L) / L
#         stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)

#         # stabilization
#         ωij = influence_function(mat, params, L) * bond_active[bond_id]

#         # update of force density
#         tij = ωij * PKinv * ΔXij
#         update_add_b_int!(storage, i, tij .* volume[j])
#         update_add_b_int!(storage, j, -tij .* volume[i])
#     end
#     return nothing
# end

function too_much_damage!(storage::RKCStorage, system::BondSystem, mat::RKCMaterial, F, i)
    if storage.damage[i] > mat.maxdmg || containsnan(F)
        # kill all bonds of this point
        storage.bond_active[each_bond_idx(system, i)] .= false
        storage.n_active_bonds[i] = 0
        return true
    end
    return false
end
