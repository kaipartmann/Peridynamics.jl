"""
    RKCRMaterial(; kernel, model, dmgmodel, monomial, epsilon, lambda, beta)

The same as the [`RKCMaterial`](@ref) but with rotation of the stress tensor for large
deformation simulations, therefore not all models are supported.

Supported models:
- `SaintVenantKirchhoff`
- `LinearElastic`

Please take a look at the [`RKCMaterial`](@ref) docs for more information about the
material, including details about the `monomial`, `epsilon`, `lambda`, and `beta`
parameters!
"""
struct RKCRMaterial{CM,K,DM,M} <: AbstractRKCMaterial{CM,NoCorrection,M}
    kernel::K
    constitutive_model::CM
    dmgmodel::DM
    epsilon::Float64
    lambda::Float64
    beta::Float64
    function RKCRMaterial(kernel::K, cm::CM, dmgmodel::DM, ::Val{M}, epsilon::Real,
                          lambda::Real, beta::Real) where {CM,K,DM,M}
        return new{CM,K,DM,M}(kernel, cm, dmgmodel, epsilon, lambda, beta)
    end
end

function RKCRMaterial(; kernel::Function=const_one_kernel,
                        model::AbstractConstitutiveModel=SaintVenantKirchhoff(),
                        dmgmodel::AbstractDamageModel=CriticalStretch(),
                        monomial::Symbol=:C1, epsilon=nothing, lambda=nothing, beta=nothing)
    if !(typeof(model) <: Union{SaintVenantKirchhoff,LinearElastic})
        msg = "model `$(typeof(model))` is currently not supported for `RKCRMaterial`!\n"
        throw(ArgumentError(msg))
    end
    get_q_dim(monomial) # check if the monomial is implemented
    ε, λ, β = get_invreg_params(epsilon, lambda, beta)
    return RKCRMaterial(kernel, model, dmgmodel, Val(monomial), ε, λ, β)
end

@params RKCRMaterial RKCPointParameters

"""
    RKCRStorage

$(extension_api_note())

Storage of [`RKCRMaterial`](@ref): the fields of [`RKCFields`](@ref) and of the fracture
bookkeeping, the rate of the deformation gradient of every point, and the left stretch, the
rotation and the unrotated Cauchy stress of every bond that the stress rotation integrates.
It carries the blocks of [`VelocityVerlet`](@ref) and [`DynamicRelaxation`](@ref) only.

$(block_table(RKCRStorage))
"""
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
    dmg_state::DamageState
end

rkc_lth_after_fields(::RKCRMaterial) = (:defgrad, :defgrad_dot, :weighted_volume)

function rkc_defgrad!(storage::RKCRStorage, system::AbstractBondSystem, mat::RKCRMaterial,
                      params::RKCPointParameters, t, Δt, i)
    (; defgrad, defgrad_dot, gradient_weight) = storage
    F = SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    Ḟ = zero(SMatrix{3,3,Float64,9})
    for bond_id in each_bond_idx(system, i)
        (; j) = get_bond(system, bond_id)
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
                              mat::RKCRMaterial, params::RKCPointParameters, t, Δt, i)
    (; volume) = system
    (; bond_active, defgrad, defgrad_dot, weighted_volume,
       bond_first_piola_kirchhoff) = storage
    Fi = get_tensor(defgrad, i)
    Ḟi = get_tensor(defgrad_dot, i)
    wi = weighted_volume[i]
    ∑P = zero(SMatrix{3,3,Float64,9})
    isolated_point(wi) && return ∑P # see `rkc_stress_integral!` of `RKCMaterial`
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            (; j, L) = get_bond(system, bond_id)
            wj = weighted_volume[j]
            if isolated_point(wj)
                update_tensor!(bond_first_piola_kirchhoff, bond_id, zero(SMatrix{3,3,Float64,9}))
                continue
            end
            ΔXij = get_vector_diff(system.position, i, j)
            Δxij = get_vector_diff(storage.position, i, j)
            Δvij = get_vector_diff(storage.velocity_half, i, j)
            Fj = get_tensor(defgrad, j)
            Ḟj = get_tensor(defgrad_dot, j)
            Fij = bond_avg(Fi, Fj, ΔXij, Δxij, L)
            Ḟij = bond_avg(Ḟi, Ḟj, ΔXij, Δvij, L)
            Pij = calc_first_piola_kirchhoff!(storage, mat, params, Fij, Ḟij, Δt, bond_id)
            Tempij = temp_ij(ΔXij, L)
            ϕ = (0.5 / wi + 0.5 / wj)
            ω̃ij = kernel(system, bond_id) * ϕ * volume[j]
            ∑Pij = ω̃ij * (Pij * Tempij)
            ∑P += ∑Pij
        end
    end
    return ∑P
end

function calc_first_piola_kirchhoff!(storage::RKCRStorage, mat::RKCRMaterial,
                                     params::RKCPointParameters, F::SMatrix{3,3,FT,9},
                                     Ḟ::SMatrix{3,3,FT,9}, Δt, bond_id) where {FT}
    D = init_stress_rotation!(storage, F, Ḟ, Δt, bond_id)
    Δε = D * Δt
    Δθ = tr(Δε)
    Δεᵈᵉᵛ = Δε - Δθ / 3 * I
    σ = get_tensor(storage.bond_unrot_cauchy_stress, bond_id)
    σₙ₊₁ = σ + 2 * params.G * Δεᵈᵉᵛ + params.K * Δθ * I
    update_tensor!(storage.bond_unrot_cauchy_stress, bond_id, σₙ₊₁)
    T = rotate_stress(storage, σₙ₊₁, bond_id)
    # the accumulated stress stays undegraded, the bond transmits the degraded stress
    P = bond_integrity(mat.dmgmodel, storage, bond_id) * first_piola_kirchhoff(T, F)
    update_tensor!(storage.bond_first_piola_kirchhoff, bond_id, P)
    return P
end

# This has to be done here, otherwise the type RKCRStorage is not known
custom_field(::Type{<:RKCRStorage}, ::Val{:hydrostatic_stress}) = true
