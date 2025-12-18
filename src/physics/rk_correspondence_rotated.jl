"""
    RKCRMaterial(; kernel, model, dmgmodel, monomial, lambda, beta)

The same as the [`RKCMaterial`](@ref) but with rotation of the stress tensor for large
deformation simulations, therefore not all models are supported.

Supported models:
- `SaintVenantKirchhoff`
- `LinearElastic`

Please take a look at the [`RKCMaterial`](@ref) docs for more information about the
material, including details about the `monomial`, `lambda`, and `beta` parameters!
"""
struct RKCRMaterial{CM,K,DM} <: AbstractRKCMaterial{CM,NoCorrection}
    kernel::K
    constitutive_model::CM
    dmgmodel::DM
    monomial::Symbol
    lambda::Float64
    beta::Float64
    function RKCRMaterial(kernel::K, cm::CM, dmgmodel::DM, monomial::Symbol,
                          lambda::Real, beta::Real) where {CM,K,DM}
        return new{CM,K,DM}(kernel, cm, dmgmodel, monomial, lambda, beta)
    end
end

function RKCRMaterial(; kernel::Function=const_one_kernel,
                        model::AbstractConstitutiveModel=SaintVenantKirchhoff(),
                        dmgmodel::AbstractDamageModel=CriticalStretch(),
                        monomial::Symbol=:C1, lambda::Real=0, beta::Real=sqrt(eps()))
    if !(typeof(model) <: Union{SaintVenantKirchhoff,LinearElastic})
        msg = "model `$(typeof(model))` is currently not supported for `RKCRMaterial`!\n"
        throw(ArgumentError(msg))
    end
    get_q_dim(monomial) # check if the kernel is implemented
    if lambda < 0
        msg = "Tikhonov regularization parameter must be non-negative! (`lambda ≥ 0`)\n"
        throw(ArgumentError(msg))
    end
    if beta < 0
        msg = "SVD truncation parameter must be non-negative! (`beta ≥ 0`)\n"
        throw(ArgumentError(msg))
    end
    return RKCRMaterial(kernel, model, dmgmodel, monomial, lambda, beta)
end

@params RKCRMaterial RKCPointParameters

@storage RKCRMaterial struct RKCRStorage
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
    @pointfield n_active_bonds::Vector{Int}
    @pointfield update_gradients::Vector{Bool}
    @pointfield cauchy_stress::Matrix{Float64}
    @pointfield von_mises_stress::Vector{Float64}
    @pointfield strain_energy_density::Vector{Float64}
    @lthfield defgrad::Matrix{Float64}
    @lthfield defgrad_dot::Matrix{Float64}
    @lthfield weighted_volume::Vector{Float64}
    bond_active::Vector{Bool}
    bond_softening::Vector{Float64}  # Gradual damage variable for stress-based softening
    gradient_weight::Matrix{Float64}
    bond_first_piola_kirchhoff::Matrix{Float64}
    left_stretch::Matrix{Float64}
    rotation::Matrix{Float64}
    bond_unrot_cauchy_stress::Matrix{Float64}
end

function init_field(::RKCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:bond_softening})
    return zeros(get_n_bonds(system))  # Initialize all bonds with zero damage
end

function init_field(::RKCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:velocity_half})
    return zeros(3, get_n_points(system))
end

function init_field(::RKCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:defgrad_dot})
    return zeros(9, get_n_points(system))
end

function init_field(::RKCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:left_stretch})
    V = zeros(9, get_n_bonds(system))
    V[[1, 5, 9], :] .= 1.0
    return V
end

function init_field(::RKCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:rotation})
    R = zeros(9, get_n_bonds(system))
    R[[1, 5, 9], :] .= 1.0
    return R
end

function init_field(::RKCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:bond_unrot_cauchy_stress})
    return zeros(9, get_n_bonds(system))
end

rkc_lth_after_fields(::RKCRMaterial) = (:defgrad, :defgrad_dot, :weighted_volume)

function rkc_defgrad!(storage::RKCRStorage, system::AbstractBondSystem, mat::RKCRMaterial,
                      params::AbstractPointParameters, t, Δt, i)
    (; bonds) = system
    (; defgrad, defgrad_dot, gradient_weight) = storage
    F = SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    Ḟ = zero(SMatrix{3,3,Float64,9})
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        Δxij = get_vector_diff(storage.position, i, j)
        Δuij = Δxij - ΔXij
        Δvij = get_vector_diff(storage.velocity_half, i, j)
        Φij = get_vector(gradient_weight, bond_id)
        F += Δuij * Φij' # maybe calculating the displacement gradient is more stable?
        Ḟ += Δvij * Φij'
    end
    update_tensor!(defgrad, i, F)
    update_tensor!(defgrad_dot, i, Ḟ)

    return nothing
end

function rkc_stress_integral!(storage::RKCRStorage, system::AbstractBondSystem,
                              mat::RKCRMaterial, params::AbstractPointParameters, t, Δt, i)
    (; bonds, volume) = system
    (; bond_active, bond_softening, defgrad, defgrad_dot, weighted_volume) = storage
    Fi = get_tensor(defgrad, i)
    Ḟi = get_tensor(defgrad_dot, i)
    wi = weighted_volume[i]

    # Early return if weighted volume is too small (fully fragmented point)
    if wi < 1e-30
        return zero(SMatrix{3,3,Float64,9})
    end
    # Get parameters from damage model
    use_tensile_split = _use_tensile_split(mat.dmgmodel)
    g_d_min = _get_residual_stiffness(mat.dmgmodel)
    ∑P = zero(SMatrix{3,3,Float64,9})
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length
            ΔXij = get_vector_diff(system.position, i, j)
            Δxij = get_vector_diff(storage.position, i, j)
            Δvij = get_vector_diff(storage.velocity_half, i, j)
            Fj = get_tensor(defgrad, j)
            Ḟj = get_tensor(defgrad_dot, j)
            Fij = bond_avg(Fi, Fj, ΔXij, Δxij, L)
            Ḟij = bond_avg(Ḟi, Ḟj, ΔXij, Δvij, L)
            Pij = calc_first_piola_kirchhoff!(storage, mat, params, Fij, Ḟij, Δt, bond_id)
            # Skip if stress calculation produced NaN
            if any(isnan, Pij) || any(isinf, Pij)
                continue
            end
                        Tempij = I - (ΔXij * ΔXij') / (L * L)
            wj = max(weighted_volume[j], 1e-30)  # Prevent division by zero
            ϕ = (0.5 / wi + 0.5 / wj)

            d = bond_softening[bond_id]
            ωij_base = kernel(system, bond_id) * ϕ * volume[j]

            if use_tensile_split && d > 0.0
                # Phase-field style: degrade only tensile stress, keep compressive intact
                Pij_eff = degraded_stress(Pij, Fij, d)
                ∑Pij = ωij_base * (Pij_eff * Tempij)
            else
                # Standard degradation: g(d) = max(ε, (1-d)²) with stiffness floor
                degradation = max(g_d_min, (1.0 - d) * (1.0 - d))
                ω̃ij = ωij_base * degradation
                ∑Pij = ω̃ij * (Pij * Tempij)
            end

            # Safety check: skip if result has NaN
            if !any(isnan, ∑Pij) && !any(isinf, ∑Pij)
                ∑P += ∑Pij
            end
        end
    end
    return ∑P
end

function calc_first_piola_kirchhoff!(storage::RKCRStorage, mat::RKCRMaterial,
                                     params::AbstractPointParameters, F::SMatrix{3,3,FT,9},
                                     Ḟ::SMatrix{3,3,FT,9}, Δt, bond_id) where {FT}
    D = init_stress_rotation!(storage, F, Ḟ, Δt, bond_id)
    Δε = D * Δt
    Δθ = tr(Δε)
    Δεᵈᵉᵛ = Δε - Δθ / 3 * I
    σ = get_tensor(storage.bond_unrot_cauchy_stress, bond_id)
    σₙ₊₁ = σ + 2 * params.G * Δεᵈᵉᵛ + params.K * Δθ * I
    update_tensor!(storage.bond_unrot_cauchy_stress, bond_id, σₙ₊₁)
    T = rotate_stress(storage, σₙ₊₁, bond_id)
    P = first_piola_kirchhoff(T, F)
    update_tensor!(storage.bond_first_piola_kirchhoff, bond_id, P)
    return P
end

# This has to be done here, otherwise the type RKCRStorage is not known
custom_field(::Type{RKCRStorage}, ::Val{:hydrostatic_stress}) = true

#-------------------------------------------------------------------------------------------
# EquivalentStress damage implementation for RKCRMaterial
#-------------------------------------------------------------------------------------------

"""
    calc_failure!(storage::RKCRStorage, system, mat, dmgmodel::EquivalentStress, paramsetup, i)

Calculate bond failure for the EquivalentStress damage model with RKCRMaterial.
Uses the stress from the previous timestep (stored in `bond_first_piola_kirchhoff`).
"""
function calc_failure!(storage::RKCRStorage, system::AbstractBondSystem,
                       mat::RKCRMaterial, dmgmodel::EquivalentStress,
                       paramsetup::AbstractParameterSetup, i)
    params = get_params(paramsetup, i)
    σc = params.σc
    (; position, n_active_bonds, bond_active, bond_softening,
       bond_first_piola_kirchhoff, defgrad) = storage
    (; bonds) = system

    # Skip if no fracture parameters defined
    if isapprox(σc, 0; atol=eps())
        for bond_id in each_bond_idx(system, i)
            n_active_bonds[i] += bond_active[bond_id]
        end
        return nothing
    end

    Fi = get_tensor(defgrad, i)

    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length

        if bond_active[bond_id]
            # Get bond deformation gradient (averaged)
            ΔXij = get_vector_diff(system.position, i, j)
            Δxij = get_vector_diff(position, i, j)
            Fj = get_tensor(defgrad, j)
            Fij = bond_avg(Fi, Fj, ΔXij, Δxij, L)

            # Get bond first Piola-Kirchhoff stress from PREVIOUS timestep and convert to Cauchy
            Pij = get_tensor(bond_first_piola_kirchhoff, bond_id)
            σij = cauchy_stress_safe(Pij, Fij)

            # Skip if stress calculation failed (NaN protection)
            if any(isnan, σij)
                n_active_bonds[i] += bond_active[bond_id]
                continue
            end

            # Calculate equivalent stress using the specified stress function
            σeq = dmgmodel.stress_func(σij)

            # Skip only if equivalent stress is NaN (numerical failure)
            if isnan(σeq)
                n_active_bonds[i] += bond_active[bond_id]
                continue
            end

            if σeq > σc && bond.fail_permit
                if dmgmodel.softening
                    # Gradual softening with rate limiting
                    softening_rate = 0.1
                    Δd_proposed = softening_rate * (σeq - σc) / σc
                    Δd = min(Δd_proposed, dmgmodel.max_damage_rate)
                    bond_softening[bond_id] = min(1.0, bond_softening[bond_id] + Δd)

                    if bond_softening[bond_id] >= 1.0
                        bond_active[bond_id] = false
                        storage.update_gradients[i] = true
                    end
                else
                    # Immediate failure
                    bond_active[bond_id] = false
                    bond_softening[bond_id] = 1.0
                    storage.update_gradients[i] = true
                end
            end
        end

        n_active_bonds[i] += bond_active[bond_id]
    end

    return nothing
end

"""
    calc_damage!(storage::RKCRStorage, system, mat, dmgmodel::EquivalentStress, paramsetup, i)

Calculate point damage for EquivalentStress with RKCRMaterial.
"""
function calc_damage!(storage::RKCRStorage, system::AbstractBondSystem,
                      mat::RKCRMaterial, dmgmodel::EquivalentStress,
                      paramsetup::AbstractParameterSetup, i)
    (; n_neighbors) = system
    (; n_active_bonds, damage, update_gradients, bond_softening, bond_active) = storage

    total_softening = 0.0
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            total_softening += bond_softening[bond_id]
        else
            total_softening += 1.0
        end
    end

    old_damage = damage[i]
    new_damage = total_softening / n_neighbors[i]

    if new_damage > old_damage + 1e-10
        update_gradients[i] = true
    end

    damage[i] = new_damage
    return nothing
end
