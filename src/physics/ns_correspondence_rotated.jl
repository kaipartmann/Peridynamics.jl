"""
    NSCRMaterial(; maxdmg, zem)

A material type used to assign the material of a [`Body`](@ref) with the local continuum
consistent (correspondence) formulation of non-ordinary state-based peridynamics.

# Keywords
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. If this value is
    exceeded, all bonds of that point are broken because the deformation gradient would then
    possibly contain `NaN` values.
    (default: `0.95`)

# Examples

```julia-repl
julia> mat = NSCRMaterial()
NSCRMaterial(maxdmg=0.95)
```

---

```julia
NSCRMaterial
```

Material type for the local continuum consistent (correspondence) formulation of
non-ordinary state-based peridynamics.

# Fields
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. See the
    constructor docs for more informations.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `NSCRMaterial`, then the following
parameters are allowed:
- `horizon::Float64`: Radius of point interactions
- `rho::Float64`: Density
- `E::Float64`: Young's modulus
- `nu::Float64`: Poisson's ratio
- `Gc::Float64`: Critical energy release rate
- `epsilon_c::Float64`: Critical strain

# Allowed export fields
When specifying the `fields` keyword of [`Job`](@ref) for a [`Body`](@ref) with
`NSCRMaterial`, the following fields are allowed:
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
struct NSCRMaterial{CM,K} <: AbstractCorrespondenceMaterial{CM,NoCorrection}
    kernel::K
    constitutive_model::CM
    maxdmg::Float64
    function NSCRMaterial(kernel::K, cm::CM, maxdmg::Real) where {CM,K}
        return new{CM,K}(kernel, cm, maxdmg)
    end
end

function Base.show(io::IO, @nospecialize(mat::NSCRMaterial))
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat, (:maxdmg,)))
    return nothing
end

function NSCRMaterial(; kernel::Function=cubic_b_spline,
                      model::AbstractConstitutiveModel=SaintVenantKirchhoff(),
                      maxdmg::Real=0.85)
    return NSCRMaterial(kernel, model, maxdmg)
end

struct NSCRPointParameters <: AbstractPointParameters
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

function NSCRPointParameters(mat::NSCRMaterial, p::Dict{Symbol,Any})
    (; δ, rho, E, nu, G, K, λ, μ) = get_required_point_parameters(mat, p)
    (; Gc, εc) = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return NSCRPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params NSCRMaterial NSCRPointParameters

@storage NSCRMaterial struct NSCRStorage
    @lthfield position::Matrix{Float64}
    @pointfield displacement::Matrix{Float64}
    @pointfield velocity::Matrix{Float64}
    @lthfield velocity_half::Matrix{Float64}
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
    # @pointfield stress::Matrix{Float64}
    # @pointfield von_mises_stress::Vector{Float64}
    @lthfield defgrad::Matrix{Float64}
    @lthfield defgrad_dot::Matrix{Float64}
    @lthfield weighted_volume::Vector{Float64}
    gradient_weight::Matrix{Float64}
    rotation::Matrix{Float64}
    left_stretch::Matrix{Float64}
    unrot_rate_of_deformation::Matrix{Float64}
    bond_stress::Matrix{Float64}
    first_piola_kirchhoff::Matrix{Float64}
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:velocity_half})
    return zeros(3, get_n_points(system))
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem, ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:damage_changed})
    return ones(Bool, get_n_loc_points(system))
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:stress})
    return zeros(9, get_n_loc_points(system))
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:von_mises_stress})
    return zeros(get_n_loc_points(system))
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:rotation})
    R = zeros(9, get_n_bonds(system))
    R[[1, 5, 9], :] .= 1.0
    return R
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:left_stretch})
    V = zeros(9, get_n_bonds(system))
    V[[1, 5, 9], :] .= 1.0
    return V
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:unrot_rate_of_deformation})
    return zeros(9, get_n_bonds(system))
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:defgrad})
    return zeros(9, get_n_points(system))
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:defgrad_dot})
    return zeros(9, get_n_points(system))
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:weighted_volume})
    return zeros(get_n_points(system))
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:gradient_weight})
    return zeros(3, get_n_bonds(system))
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:first_piola_kirchhoff})
    return zeros(9, get_n_bonds(system))
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:bond_stress})
    return zeros(9, get_n_bonds(system))
end

function initialize!(chunk::BodyChunk{<:BondSystem,<:NSCRMaterial})
    chunk.storage.damage_changed .= true
    return nothing
end

@inline function calc_damage!(chunk::BodyChunk{<:BondSystem,<:NSCRMaterial})
    (; system, storage, paramsetup) = chunk
    (; n_neighbors, bonds) = system
    (; n_active_bonds, damage, damage_changed) = storage
    for i in each_point_idx(chunk)
        params = get_params(paramsetup, i)
        for bond_id in each_bond_idx(system, i)
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length
            Δxij = get_coordinates_diff(storage, i, j)
            l = norm(Δxij)
            ε = (l - L) / L
            stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)
        end
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

function calc_force_density!(dh::ThreadsBodyDataHandler{<:BondSystem,<:NSCRMaterial}, t, Δt)
    @threads :static for chunk_id in eachindex(dh.chunks)
        exchange_loc_to_halo!(dh, chunk_id, (:position, :velocity_half))
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        calc_weights_and_defgrad!(dh.chunks[chunk_id], t, Δt)
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        exchange_loc_to_halo!(dh, chunk_id, (:defgrad, :defgrad_dot, :weighted_volume))
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

function calc_force_density!(dh::MPIBodyDataHandler{<:BondSystem,<:NSCRMaterial}, t, Δt)
    (; chunk) = dh
    exchange_loc_to_halo!(dh, (:position, :velocity_half))
    calc_weights_and_defgrad!(chunk, t, Δt)
    exchange_loc_to_halo!(dh, (:defgrad, :defgrad_dot, :weighted_volume))
    calc_force_density!(chunk, t, Δt)
    nancheck(chunk, t)
    exchange_halo_to_loc!(dh)
    return nothing
end

function calc_weights_and_defgrad!(chunk::BodyChunk{<:BondSystem,<:NSCRMaterial}, t, Δt)
    (; system, mat, paramsetup, storage) = chunk
    for i in each_point_idx(system)
        calc_failure_point!(storage, system, mat, paramsetup, i)
        calc_weights_and_defgrad!(storage, system, mat, paramsetup, t, Δt, i)
    end
    return nothing
end

function calc_weights_and_defgrad!(storage::NSCRStorage, system::BondSystem,
                                   mat::NSCRMaterial,
                                   paramhandler::AbstractParameterHandler, t, Δt, i)
    params = get_params(paramhandler, i)
    calc_weights_and_defgrad!(storage, system, mat, params, t, Δt, i)
    return nothing
end

function calc_weights_and_defgrad!(storage::NSCRStorage, system::BondSystem,
                                   mat::NSCRMaterial, params::NSCRPointParameters, t, Δt, i)
    (; bonds, volume) = system
    (; bond_active, weighted_volume, gradient_weight, damage_changed) = storage

    K = zero(SMatrix{3,3,Float64,9})
    _F = zero(SMatrix{3,3,Float64,9})
    _Ḟ = zero(SMatrix{3,3,Float64,9})
    wi = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_diff(system.position, i, j)
        Δxij = get_diff(storage.position, i, j)
        Δvij = get_diff(storage.velocity_half, i, j)
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        temp = ωij * volume[j]
        K += temp * (ΔXij * ΔXij')
        _F += temp * (Δxij * ΔXij')
        _Ḟ += temp * (Δvij * ΔXij')
        wi += temp
    end
    Kinv = inv(K)
    F = _F * Kinv
    Ḟ = _Ḟ * Kinv
    update_tensor!(storage.defgrad, i, F)
    update_tensor!(storage.defgrad_dot, i, Ḟ)
    weighted_volume[i] = wi

    damage_changed[i] || return nothing

    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_diff(system.position, i, j)
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        temp = ωij * volume[j]
        Ψ = temp * (Kinv * ΔXij)
        update_vector!(gradient_weight, bond_id, Ψ)
    end

    return nothing
end

function force_density_point!(storage::NSCRStorage, system::BondSystem, mat::NSCRMaterial,
                              paramhandler::AbstractParameterHandler, t, Δt, i)
    params = get_params(paramhandler, i)
    force_density_point!(storage, system, mat, params, t, Δt, i)
    return nothing
end

function force_density_point!(storage::NSCRStorage, system::BondSystem, mat::NSCRMaterial,
                              params::NSCRPointParameters, t, Δt, i)
    # Fi = get_tensor(storage.defgrad, i)
    # too_much_damage!(storage, system, mat, Fi, i) && return nothing
    # P = calc_first_piola_kirchhoff!(storage, mat, params, Fi, Δt, i)
    # nsc_force_density!(storage, system, mat, params, P, i)
    # nsc_force_density!(storage, system, mat, params, i)

    (; bonds, volume) = system
    (; bond_active, gradient_weight, defgrad, defgrad_dot, weighted_volume) = storage

    Fi = get_tensor(defgrad, i)
    Ḟi = get_tensor(defgrad_dot, i)
    # too_much_damage!(storage, system, mat, Fi, i) && return nothing

    # Stress integral ∑P
    wi = weighted_volume[i]
    ∑P = SMatrix{3,3,Float64,9}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_coordinates_diff(storage, i, j)
        # l = norm(Δxij)
        # ε = (l - L) / L
        # stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)

        if bond_active[bond_id]
            ΔXij = get_coordinates_diff(system, i, j)
            Δvij = get_diff(storage.velocity_half, i, j)
            Fj = get_tensor(defgrad, j)
            Ḟj = get_tensor(defgrad_dot, j)
            Fb = 0.5 * (Fi + Fj)
            Ḟb = 0.5 * (Ḟi + Ḟj)
            ΔXijLL = ΔXij' / (L * L)
            ΔFij = (Δxij - Fb * ΔXij) * ΔXijLL
            ΔḞij = (Δvij - Ḟb * ΔXij) * ΔXijLL
            Fij = Fb + ΔFij
            Ḟij = Ḟb + ΔḞij
            Pij = calc_first_piola_kirchhoff!(storage, mat, params, Fij, Ḟij, Δt, bond_id)
            Tempij = I - ΔXij * ΔXijLL
            # ω̃ij = kernel(system, bond_id) * (0.5 / wi + 0.5 / weighted_volume[j])
            wj = weighted_volume[j]
            ϕ = (wi > 0 && wj > 0) ? (0.5 / wi + 0.5 / wj) : 0.0
            ω̃ij = kernel(system, bond_id) * ϕ
            ∑P += ω̃ij * (Pij * Tempij)
        end
    end

    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length
            ΔXij = get_coordinates_diff(system, i, j)
            Pij = get_tensor(storage.first_piola_kirchhoff, bond_id)
            Φij = get_vector(gradient_weight, bond_id)
            # ω̃ij = kernel(system, bond_id) * (0.5 / wi + 0.5 / weighted_volume[j])
            wj = weighted_volume[j]
            ϕ = (wi > 0 && wj > 0) ? (0.5 / wi + 0.5 / wj) : 0.0
            ω̃ij = kernel(system, bond_id) * ϕ
            tij = ω̃ij / (L * L) * (Pij * ΔXij) + ∑P * Φij
            update_add_b_int!(storage, i, tij * volume[j])
            update_add_b_int!(storage, j, -tij * volume[i])
        end
    end

    return nothing
end

# function calc_first_piola_kirchhoff!(storage::NSCRStorage, mat::NSCRMaterial,
#                                      params::NSCRPointParameters, F, Ḟ, Δt, i)
#     P = first_piola_kirchhoff(mat.constitutive_model, storage, params, F)
#     σ = cauchy_stress(P, F)
#     init_stress_rotation!(storage, F, Ḟ, Δt, i)
#     T = rotate_stress(storage, σ, i)
#     P = first_piola_kirchhoff(T, F)
#     update_tensor!(storage.first_piola_kirchhoff, i, P)
#     return P
# end

function calc_first_piola_kirchhoff!(storage::NSCRStorage, ::NSCRMaterial,
                                     params::NSCRPointParameters, F, Ḟ, Δt, bond_id)
    D = init_stress_rotation!(storage, F, Ḟ, Δt, bond_id)
    Δε = D * Δt
    Δθ = tr(Δε)
    Δεᵈᵉᵛ = Δε - Δθ / 3 * I
    σ = get_tensor(storage.bond_stress, bond_id)
    σₙ₊₁ = σ + 2 * params.G * Δεᵈᵉᵛ + params.K * Δθ * I
    update_tensor!(storage.bond_stress, bond_id, σₙ₊₁)
    T = rotate_stress(storage, σₙ₊₁, bond_id)
    P = first_piola_kirchhoff(T, F)
    update_tensor!(storage.first_piola_kirchhoff, bond_id, P)
    return P
end

# function calc_first_piola_kirchhoff!(storage, mat, params, F, Ḟ, Δt, i)
#     𝐋 = @MMatrix zeros(3,3)
#     RoD = @MMatrix zeros(3,3)
#     Spin = @MMatrix zeros(3,3)
#     left_stretch = @MMatrix zeros(3,3)
#     z = @MVector zeros(3)
#     w = @MVector zeros(3)
#     𝛀Tens = @MMatrix zeros(3,3)
#     Qmatrix = @MMatrix zeros(3,3)
#     rot_tens_ = @MMatrix zeros(3,3)
#     𝛔_unrot = @MMatrix zeros(3,3)
#     tempVec = @MMatrix zeros(3,3)
#     𝐅 = F
#     𝐅dot = Ḟ
#     𝐋 = 𝐅dot * inv(𝐅)
#     RoD = 1/2 * (𝐋 + transpose(𝐋))
#     Spin = 1/2 * (𝐋 - transpose(𝐋))
#     left_stretch[1,1] = storage.left_stretch[1,i]
#     left_stretch[1,2] = storage.left_stretch[2,i]
#     left_stretch[1,3] = storage.left_stretch[3,i]
#     left_stretch[2,1] = storage.left_stretch[4,i]
#     left_stretch[2,2] = storage.left_stretch[5,i]
#     left_stretch[2,3] = storage.left_stretch[6,i]
#     left_stretch[3,1] = storage.left_stretch[7,i]
#     left_stretch[3,2] = storage.left_stretch[8,i]
#     left_stretch[3,3] = storage.left_stretch[9,i]
#     z[1] = - storage.left_stretch[3, i] * RoD[4] - storage.left_stretch[6, i] * RoD[5] -
#              storage.left_stretch[9, i] * RoD[6] + storage.left_stretch[2, i] * RoD[7] +
#              storage.left_stretch[5, i] * RoD[8] + storage.left_stretch[8, i] * RoD[9];
#     z[2] =   storage.left_stretch[3, i] * RoD[1] + storage.left_stretch[6, i] * RoD[2] +
#              storage.left_stretch[9, i] * RoD[3] - storage.left_stretch[1, i] * RoD[7] -
#              storage.left_stretch[4, i] * RoD[8] - storage.left_stretch[7, i] * RoD[9];
#     z[3] = - storage.left_stretch[2, i] * RoD[1] - storage.left_stretch[5, i] * RoD[2] -
#              storage.left_stretch[8, i] * RoD[3] + storage.left_stretch[1, i] * RoD[4] +
#              storage.left_stretch[4, i] * RoD[5] + storage.left_stretch[7, i] * RoD[6];
#     w[1] = 0.5 * (Spin[3,2] - Spin[2,3])
#     w[2] = 0.5 * (Spin[1,3] - Spin[3,1])
#     w[3] = 0.5 * (Spin[2,1] - Spin[1,2])
#     traceV = storage.left_stretch[1, i] + storage.left_stretch[5, i] + storage.left_stretch[9, i]
#     omega = w + inv(traceV * I - left_stretch) * z
#     𝛀Tens[1,1] = 0.0
#     𝛀Tens[1,2] = -omega[3]
#     𝛀Tens[1,3] = omega[2]
#     𝛀Tens[2,1] = omega[3]
#     𝛀Tens[2,2] = 0.0
#     𝛀Tens[2,3] = -omega[1]
#     𝛀Tens[3,1] = -omega[2]
#     𝛀Tens[3,2] = omega[1]
#     𝛀Tens[3,3] = 0.0
#     ΩSq = omega[1]^2 + omega[2]^2 + omega[3]^2
#     Ω = sqrt(ΩSq)
#     if ΩSq > 1e-16 && Ω !== Inf
#         scfac1 = sin(Δt * Ω) / Ω
#         scfac2 = -(1 - cos(Δt * Ω)) / ΩSq
#         𝛀TensSq = 𝛀Tens * 𝛀Tens
#         Qmatrix = I + scfac1 * 𝛀Tens + scfac2 * 𝛀TensSq
#     else
#         Qmatrix .= 0.
#         Qmatrix[1,1] = 1.0
#         Qmatrix[2,2] = 1.0
#         Qmatrix[3,3] = 1.0
#     end
#     rot_tens_[1,1] = storage.rotation[1,i]
#     rot_tens_[1,2] = storage.rotation[2,i]
#     rot_tens_[1,3] = storage.rotation[3,i]
#     rot_tens_[2,1] = storage.rotation[4,i]
#     rot_tens_[2,2] = storage.rotation[5,i]
#     rot_tens_[2,3] = storage.rotation[6,i]
#     rot_tens_[3,1] = storage.rotation[7,i]
#     rot_tens_[3,2] = storage.rotation[8,i]
#     rot_tens_[3,3] = storage.rotation[9,i]
#     rot_tens = Qmatrix * rot_tens_
#     storage.rotation[1,i] = rot_tens[1,1]
#     storage.rotation[2,i] = rot_tens[1,2]
#     storage.rotation[3,i] = rot_tens[1,3]
#     storage.rotation[4,i] = rot_tens[2,1]
#     storage.rotation[5,i] = rot_tens[2,2]
#     storage.rotation[6,i] = rot_tens[2,3]
#     storage.rotation[7,i] = rot_tens[3,1]
#     storage.rotation[8,i] = rot_tens[3,2]
#     storage.rotation[9,i] = rot_tens[3,3]
#     Vdot = 𝐋 * left_stretch - left_stretch * 𝛀Tens
#     storage.left_stretch[1,i] += Δt * Vdot[1,1]
#     storage.left_stretch[2,i] += Δt * Vdot[1,2]
#     storage.left_stretch[3,i] += Δt * Vdot[1,3]
#     storage.left_stretch[4,i] += Δt * Vdot[2,1]
#     storage.left_stretch[5,i] += Δt * Vdot[2,2]
#     storage.left_stretch[6,i] += Δt * Vdot[2,3]
#     storage.left_stretch[7,i] += Δt * Vdot[3,1]
#     storage.left_stretch[8,i] += Δt * Vdot[3,2]
#     storage.left_stretch[9,i] += Δt * Vdot[3,3]
#     tempVec = RoD * rot_tens
#     UnRotRoD = transpose(rot_tens) * tempVec
#     strainInc = UnRotRoD * Δt
#     deviatoricStrain = copy(strainInc)
#     dilatation = strainInc[1,1] + strainInc[2,2] + strainInc[3,3]
#     deviatoricStrain[1,1] -= dilatation/3
#     deviatoricStrain[2,2] -= dilatation/3
#     deviatoricStrain[3,3] -= dilatation/3
#     storage.bond_stress[1,i] += deviatoricStrain[1,1] * 2 * params.G + params.K * dilatation
#     storage.bond_stress[2,i] += deviatoricStrain[1,2] * 2 * params.G
#     storage.bond_stress[3,i] += deviatoricStrain[1,3] * 2 * params.G
#     storage.bond_stress[4,i] += deviatoricStrain[2,1] * 2 * params.G
#     storage.bond_stress[5,i] += deviatoricStrain[2,2] * 2 * params.G + params.K * dilatation
#     storage.bond_stress[6,i] += deviatoricStrain[2,3] * 2 * params.G
#     storage.bond_stress[7,i] += deviatoricStrain[3,1] * 2 * params.G
#     storage.bond_stress[8,i] += deviatoricStrain[3,2] * 2 * params.G
#     storage.bond_stress[9,i] += deviatoricStrain[3,3] * 2 * params.G + params.K * dilatation
#     𝛔_unrot[1,1] = storage.bond_stress[1,i]
#     𝛔_unrot[1,2] = storage.bond_stress[2,i]
#     𝛔_unrot[1,3] = storage.bond_stress[3,i]
#     𝛔_unrot[2,1] = storage.bond_stress[4,i]
#     𝛔_unrot[2,2] = storage.bond_stress[5,i]
#     𝛔_unrot[2,3] = storage.bond_stress[6,i]
#     𝛔_unrot[3,1] = storage.bond_stress[7,i]
#     𝛔_unrot[3,2] = storage.bond_stress[8,i]
#     𝛔_unrot[3,3] = storage.bond_stress[9,i]
#     tempVec = 𝛔_unrot * transpose(rot_tens)
#     T = rot_tens * tempVec
#     J = det(F)
#     P = J * T * inv(F)'
#     update_tensor!(storage.first_piola_kirchhoff, i, SMatrix{3,3,Float64,9}(P))
#     return P
# end

# Working version!!!
# function nsc_force_density!(storage::NSCRStorage, system::BondSystem, mat::NSCRMaterial,
#                             params::NSCRPointParameters, P::SMatrix, i)
#     (; bonds, volume) = system
#     (; bond_active, gradient_weight) = storage
#     for bond_id in each_bond_idx(system, i)
#         bond = bonds[bond_id]
#         j, L = bond.neighbor, bond.length
#         Δxij = get_coordinates_diff(storage, i, j)
#         l = norm(Δxij)
#         ε = (l - L) / L
#         stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)
#         if bond_active[bond_id]
#             Ψ = get_vector(gradient_weight, bond_id)
#             tij = P * Ψ
#             update_add_b_int!(storage, i, tij .* volume[j])
#             update_add_b_int!(storage, j, -tij .* volume[i])
#         end
#     end
#     return nothing
# end

# function nsc_force_density!(storage::NSCRStorage, system::BondSystem, mat::NSCRMaterial,
#                             params::NSCRPointParameters, i)
#     (; bonds, volume) = system
#     (; bond_active, gradient_weight, defgrad) = storage
#     Fi = get_tensor(defgrad, i)
#     too_much_damage!(storage, system, mat, Fi, i) && return nothing
#     for bond_id in each_bond_idx(system, i)
#         bond = bonds[bond_id]
#         j, L = bond.neighbor, bond.length
#         Δxij = get_coordinates_diff(storage, i, j)
#         l = norm(Δxij)
#         ε = (l - L) / L
#         stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)
#         if bond_active[bond_id]
#             Fj = get_tensor(defgrad, j)
#             Fb = 0.5 * (Fi + Fj)
#             ΔXij = get_coordinates_diff(system, i, j)
#             ΔFij = (Δxij - Fb * ΔXij) * ΔXij' / (L * L)
#             Fij = Fb + ΔFij
#             P = calc_first_piola_kirchhoff(storage, mat, params, Fij)
#             Ψ = get_vector(gradient_weight, bond_id)
#             tij = P * Ψ
#             update_add_b_int!(storage, i, tij .* volume[j])
#             update_add_b_int!(storage, j, -tij .* volume[i])
#         end
#     end
#     return nothing
# end

# function nsc_force_density!(storage::NSCRStorage, system::BondSystem, mat::NSCRMaterial,
#                             params::NSCRPointParameters, i)
#     (; bonds, volume) = system
#     (; bond_active, gradient_weight, defgrad, weighted_volume) = storage

#     Fi = get_tensor(defgrad, i)
#     too_much_damage!(storage, system, mat, Fi, i) && return nothing

#     # Stress integral SI
#     wi = weighted_volume[i]
#     SI = SMatrix{3,3,Float64,9}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
#     for bond_id in each_bond_idx(system, i)
#         bond = bonds[bond_id]
#         j, L = bond.neighbor, bond.length
#         ΔXij = get_coordinates_diff(system, i, j)
#         Δxij = get_coordinates_diff(storage, i, j)
#         l = norm(Δxij)
#         ε = (l - L) / L
#         stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)

#         if bond_active[bond_id]
#             Fj = get_tensor(defgrad, j)
#             Fb = 0.5 * (Fi + Fj)
#             ΔFij = (Δxij - Fb * ΔXij) * ΔXij' / (L * L)
#             Fij = Fb + ΔFij
#             Pij = calc_first_piola_kirchhoff(storage, mat, params, Fij)
#             update_tensor!(storage.first_piola_kirchhoff, bond_id, Pij)
#             Tempij = I - (ΔXij * ΔXij') / (L * L)
#             ωij = kernel(system, bond_id)
#             wj = weighted_volume[j]
#             ϕ = volume[j] * ωij * (0.5 / wi + 0.5 / wj)
#             SI += ϕ * (Pij * Tempij)
#         end
#     end

#     for bond_id in each_bond_idx(system, i)
#         if bond_active[bond_id]
#             bond = bonds[bond_id]
#             j, L = bond.neighbor, bond.length
#             ΔXij = get_coordinates_diff(system, i, j)
#             Pij = get_tensor(storage.first_piola_kirchhoff, bond_id)
#             Φij = get_vector(gradient_weight, bond_id)
#             ωij = kernel(system, bond_id)
#             tij = ωij / (wi * L * L) * (Pij * ΔXij) + SI * Φij
#             update_add_b_int!(storage, i, tij * volume[j])
#             update_add_b_int!(storage, j, -tij * volume[i])
#         end
#     end
#     return nothing
# end

# NODALLY STABLIZIED VERSION!!! WORKING, BUT NOT SURE IF CORRECT!
# function nsc_force_density!(storage::NSCRStorage, system::BondSystem, mat::NSCRMaterial,
#                             params::NSCRPointParameters, P::SMatrix, i)
#     (; bonds, volume) = system
#     (; bond_active, gradient_weight, defgrad) = storage
#     Fi = get_tensor(defgrad, i)
#     too_much_damage!(storage, system, mat, Fi, i) && return nothing
#     for bond_id in each_bond_idx(system, i)
#         bond = bonds[bond_id]
#         j, L = bond.neighbor, bond.length
#         Δxij = get_coordinates_diff(storage, i, j)
#         l = norm(Δxij)
#         ε = (l - L) / L
#         stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)
#         if bond_active[bond_id]
#             ΔXij = get_coordinates_diff(system, i, j)
#             Fj = get_tensor(defgrad, j)
#             Fb = 0.5 * (Fi + Fj)
#             ΔFij = (Δxij - Fb * ΔXij) * ΔXij' / (L * L)
#             Fij = Fb + ΔFij
#             Pij = calc_first_piola_kirchhoff(storage, mat, params, Fij)
#             Ψ = get_vector(gradient_weight, bond_id)
#             tij = (Pij * Ψ) / volume[j]
#             update_add_b_int!(storage, i, tij .* volume[j])
#             update_add_b_int!(storage, j, -tij .* volume[i])
#         end
#     end
#     return nothing
# end


# PERIDIGM VERSION
# function nsc_force_density!(storage::NSCRStorage, system::BondSystem, mat::NSCRMaterial,
#                             params::NSCRPointParameters, P::SMatrix, i)
#     (; bonds, volume) = system
#     (; bond_active, gradient_weight, defgrad, weighted_volume) = storage
#     Fi = get_tensor(defgrad, i)
#     wi = weighted_volume[i]
#     SI = SMatrix{3,3,Float64,9}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
#     for bond_id in each_bond_idx(system, i)
#         bond = bonds[bond_id]
#         j, L = bond.neighbor, bond.length
#         ΔXij = get_coordinates_diff(system, i, j)
#         Δxij = get_coordinates_diff(storage, i, j)
#         l = norm(Δxij)
#         ε = (l - L) / L
#         stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)

#         Fj = get_tensor(defgrad, j)
#         Fb = 0.5 * (Fi + Fj)
#         ΔFij = (Δxij - Fb * ΔXij) * ΔXij' / (L * L)
#         Fij = Fb + ΔFij
#         Pij = calc_first_piola_kirchhoff(storage, mat, params, Fij)
#         # update_tensor!(first_piola_kirchhoff, bond_id, Pij)
#         # ΔPij = Pij - Pi
#         Tempij = (I - ΔXij * ΔXij') / (L * L)
#         ωij = influence_function(mat, params, L) * bond_active[bond_id]
#         wj = weighted_volume[j]
#         ϕ = volume[j] * ωij * (0.5 / wi + 0.5 / wj)
#         SI += ϕ * (Pij * Tempij)
#     end
#     for bond_id in each_bond_idx(system, i)
#         bond = bonds[bond_id]
#         j, L = bond.neighbor, bond.length
#         ΔXij = get_coordinates_diff(system, i, j)
#         Δxij = get_coordinates_diff(storage, i, j)
#         Fj = get_tensor(defgrad, j)
#         Fb = 0.5 * (Fi + Fj)
#         ΔFij = (Δxij - Fb * ΔXij) * ΔXij' / (L * L)
#         Fij = Fb + ΔFij
#         Pij = calc_first_piola_kirchhoff(storage, mat, params, Fij)
#         Φij = get_vector(gradient_weight, bond_id)
#         ωij = influence_function(mat, params, L) * bond_active[bond_id]
#         tij = (SI * Φij) / volume[j] + ωij / (wi * L * L) * (Pij * ΔXij)
#         update_add_b_int!(storage, i, tij * volume[j])
#         update_add_b_int!(storage, j, -tij * volume[i])
#     end
#     return nothing
# end

# NODALLY STABLIZIED VERSION WRONG?
# function nsc_force_density!(storage::NSCRStorage, system::BondSystem, mat::NSCRMaterial,
#                             params::NSCRPointParameters, P::SMatrix, i)
#     (; bonds, volume) = system
#     (; bond_active, gradient_weight, defgrad, weighted_volume) = storage
#     Fi = get_tensor(defgrad, i)
#     Vi = volume[i]
#     wi = weighted_volume[i]

#     for bond_id_j in each_bond_idx(system, i)
#         bond_ij = bonds[bond_id_j]
#         j, Lij = bond_ij.neighbor, bond_ij.length
#         ΔXij = get_coordinates_diff(storage, i, j)
#         Δxij = get_coordinates_diff(storage, i, j)
#         lij = norm(Δxij)
#         ε = (lij - Lij) / Lij
#         stretch_based_failure!(storage, system, bond_ij, params, ε, i, bond_id_j)

#         ωij = influence_function(mat, params, Lij)
#         Vj, wj = volume[j], weighted_volume[j]
#         Ṽij = (Vi * Vj) / 2 * (ωij / wj + ωij / wi)
#         Ψij = get_vector(gradient_weight, bond_id_j)
#         Fj = get_tensor(defgrad, j)
#         Fij = 0.5 * (Fi + Fj)
#         Pij = calc_first_piola_kirchhoff(storage, mat, params, Fij)

#         Tcor = Ṽij * (Pij * ΔXij) / (Vj * Vi * Lij * Lij)

#         Tave = zero(SVector{3,Float64})
#         for bond_id_k in each_bond_idx(system, i)
#             bond_ik = bonds[bond_id_k]
#             k, Lik = bond_ik.neighbor, bond_ik.length
#             ωik = influence_function(mat, params, Lik)
#             Vk, wk = volume[k], weighted_volume[k]
#             Ṽik = (Vi * Vk) / 2 * (ωik / wk + ωij / wi)
#             Fk = get_tensor(defgrad, k)
#             Fik = 0.5 * (Fi + Fk)
#             Pik = calc_first_piola_kirchhoff(storage, mat, params, Fik)
#             Tave += Ṽik * (Pik * Ψij) / (Vj * Vi)
#             ΔXΔXΨ = (ΔXij * ΔXij') * Ψij
#             Tcor -= Ṽik * (Pik * ΔXΔXΨ) / (Vj * Vi * Lik * Lik)
#         end

#         T = Tave + Tcor
#         update_add_b_int!(storage, i, T .* volume[j])
#         update_add_b_int!(storage, j, -T .* volume[i])
#     end
#     return nothing
# end

    # Tcor2 = SMatrix{3,3,Float64,9}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # for bond_id in each_bond_idx(system, i)
    #     bond = bonds[bond_id]
    #     j, L = bond.neighbor, bond.length
    #     ΔXij = get_coordinates_diff(system, i, j)
    #     Δxij = get_coordinates_diff(storage, i, j)
    #     l = norm(Δxij)
    #     ε = (l - L) / L
    #     stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)

    #     if bond_active[bond_id]
    #         ωij = influence_function(mat, params, L)
    #         Vj, wj = volume[j], weighted_volume[j]
    #         Ṽij = (Vi * Vj) / 2 * (ωij / wi + ωij / wj)

    #         Fj = get_tensor(defgrad, j)
    #         Fb = 0.5 * (Fi + Fj)
    #         Pb = calc_first_piola_kirchhoff(storage, mat, params, Fb)

    #         temp = Ṽij / (Vi * Vj * L * L)
    #         Ψ = get_vector(gradient_weight, bond_id)
    #         ΔXΨ = (ΔXij * ΔXij') * Ψ
    #         Tcor2 += temp * (P * ΔXΨ)
    #     end
    # end

    # for bond_id in each_bond_idx(system, i)
    #     if bond_active[bond_id]
    #         bond = bonds[bond_id]
    #         j, L = bond.neighbor, bond.length
    #         ΔXij = get_coordinates_diff(system, i, j)
    #         Δxij = get_coordinates_diff(storage, i, j)

    #         # update of force density
    #         Ψ = get_vector(gradient_weight, bond_id)
    #         tij = (P * Ψ) / volume[j]
    #         update_add_b_int!(storage, i, tij .* volume[j])
    #         update_add_b_int!(storage, j, -tij .* volume[i])
    #     end
    # end
#     return nothing
# end

function too_much_damage!(storage::NSCRStorage, system::BondSystem, mat::NSCRMaterial, F, i)
    if storage.damage[i] > mat.maxdmg || containsnan(F)
        # kill all bonds of this point
        storage.bond_active[each_bond_idx(system, i)] .= false
        storage.n_active_bonds[i] = 0
        return true
    end
    return false
end
