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

@params RKCRMaterial StandardPointParameters

@storage RKCRMaterial struct RKCRStorage
    @inherit VelocityVerletFields DynamicRelaxationFields
    @inherit BondFracFields RKCFields
    @lth velocity_half::PointVector
    @htl b_int::PointVector
    cauchy_stress::PointTensor
    von_mises_stress::PointScalar
    strain_energy_density::PointScalar
    @lth defgrad_dot::PointTensor
    left_stretch::BondTensor = I
    rotation::BondTensor = I
    bond_unrot_cauchy_stress::BondTensor
end

rkc_lth_after_fields(::RKCRMaterial) = (:defgrad, :defgrad_dot, :weighted_volume)

function rkc_defgrad!(storage::RKCRStorage, system::AbstractBondSystem, mat::RKCRMaterial,
                      params::StandardPointParameters, t, Δt, i)
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
    σ = get_tensor(storage.bond_unrot_cauchy_stress, bond_id)
    σₙ₊₁ = σ + 2 * params.G * Δεᵈᵉᵛ + params.K * Δθ * I
    update_tensor!(storage.bond_unrot_cauchy_stress, bond_id, σₙ₊₁)
    T = rotate_stress(storage, σₙ₊₁, bond_id)
    P = first_piola_kirchhoff(T, F)
    update_tensor!(storage.bond_first_piola_kirchhoff, bond_id, P)
    return P
end

# This has to be done here, otherwise the type RKCRStorage is not known
custom_field(::Type{<:RKCRStorage}, ::Val{:hydrostatic_stress}) = true
