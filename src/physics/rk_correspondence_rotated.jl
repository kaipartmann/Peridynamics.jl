"""
    RKCRMaterial(; kernel, model, dmgmodel, maxdmg, reprkernel, regfactor)

The same as the [`RKCMaterial`](@ref) but with rotation of the stress tensor for large
deformation simulations, therefore not all models are supported.

Supported models:
- `SaintVenantKirchhoff`
- `LinearElastic`

Please take a look at the [`RKCMaterial`](@ref) docs for more information about the
material!
"""
struct RKCRMaterial{CM,K,DM} <: AbstractRKCMaterial{CM,NoCorrection}
    kernel::K
    constitutive_model::CM
    dmgmodel::DM
    maxdmg::Float64
    reprkernel::Symbol
    regfactor::Float64
    function RKCRMaterial(kernel::K, cm::CM, dmgmodel::DM, maxdmg::Real, reprkernel::Symbol,
                          regfactor::Real) where {CM,K,DM}
        return new{CM,K,DM}(kernel, cm, dmgmodel, maxdmg, reprkernel, regfactor)
    end
end

function RKCRMaterial(; kernel::Function=cubic_b_spline_kernel,
                        model::AbstractConstitutiveModel=SaintVenantKirchhoff(),
                        dmgmodel::AbstractDamageModel=CriticalStretch(), maxdmg::Real=1.0,
                        reprkernel::Symbol=:C1, regfactor::Real=1e-13)
    if !(model <: Union{SaintVenantKirchhoff,LinearElastic})
        msg = "model `$(model)` is currently not supported for `CRMaterial`!\n"
        throw(ArgumentError(msg))
    end
    get_q_dim(reprkernel) # check if the kernel is implemented
    if !(0 ≤ regfactor ≤ 1)
        msg = "Regularization factor must be in the range 0 ≤ regfactor ≤ 1\n"
        throw(ArgumentError(msg))
    end
    return RKCRMaterial(kernel, model, dmgmodel, maxdmg, reprkernel, regfactor)
end

@params RKCRMaterial StandardPointParameters

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
    bond_active::Vector{Bool}
    @pointfield n_active_bonds::Vector{Int}
    @pointfield update_gradients::Vector{Bool}
    @pointfield stress::Matrix{Float64}
    @pointfield von_mises_stress::Vector{Float64}
    @lthfield defgrad::Matrix{Float64}
    @lthfield defgrad_dot::Matrix{Float64}
    @lthfield weighted_volume::Vector{Float64}
    gradient_weight::Matrix{Float64}
    first_piola_kirchhoff::Matrix{Float64}
    left_stretch::Matrix{Float64}
    rotation::Matrix{Float64}
    bond_stress::Matrix{Float64}
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
                    ::Val{:bond_stress})
    return zeros(9, get_n_bonds(system))
end

rkc_lth_after_fields(::RKCRMaterial) = (:defgrad, :defgrad_dot, :weighted_volume)

function rkc_defgrad!(storage::RKCRStorage, system::AbstractBondSystem, mat::RKCRMaterial,
                      params::StandardPointParameters, t, Δt, i)
    (; bonds) = system
    (; defgrad, defgrad_dot, gradient_weight) = storage
    F = zero(SMatrix{3,3,Float64,9})
    Ḟ = zero(SMatrix{3,3,Float64,9})
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        Δxij = get_vector_diff(storage.position, i, j)
        Δvij = get_vector_diff(storage.velocity_half, i, j)
        Φij = get_vector(gradient_weight, bond_id)
        F += Δxij * Φij'
        Ḟ += Δvij * Φij'
    end
    update_tensor!(defgrad, i, F)
    update_tensor!(defgrad_dot, i, Ḟ)

    return nothing
end

function rkc_stress_integral!(storage::RKCRStorage, system::AbstractBondSystem,
                              mat::RKCRMaterial, params::StandardPointParameters, t, Δt, i)
    (; bonds, volume) = system
    (; bond_active, defgrad, defgrad_dot, weighted_volume) = storage
    Fi = get_tensor(defgrad, i)
    Ḟi = get_tensor(defgrad_dot, i)
    wi = weighted_volume[i]
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
            Tempij = I - (ΔXij * ΔXij') / (L * L)
            wj = weighted_volume[j]
            ϕ = (0.5 / wi + 0.5 / wj)
            ω̃ij = kernel(system, bond_id) * ϕ * volume[j]
            ∑Pij = ω̃ij * (Pij * Tempij)
            ∑P += ∑Pij
        end
    end
    return ∑P
end

function calc_first_piola_kirchhoff!(storage::RKCRStorage, mat::RKCRMaterial,
                                     params::StandardPointParameters, F::SMatrix{3,3,FT,9},
                                     Ḟ::SMatrix{3,3,FT,9}, Δt, bond_id) where {FT}
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
