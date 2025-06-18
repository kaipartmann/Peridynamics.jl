"""
    RKCRMaterial(; kernel, model, dmgmodel, maxdmg, regfactor)

The same as the [`RKCMaterial`](@ref) but with rotation of the stress tensor for large
deformation simulations, therefore not all models are supported.

Supported models:
- `LinearElastic`

Please take a look at the [`RKCMaterial`](@ref) docs for more information about the
material!
"""
struct RKCRMaterial{CM,K,DM} <: AbstractRKCMaterial{CM,NoCorrection}
    kernel::K
    constitutive_model::CM
    dmgmodel::DM
    maxdmg::Float64
    accuracy_order::Int
    regfactor::Float64
    function RKCRMaterial(kernel::K, cm::CM, dmgmodel::DM, maxdmg::Real,
                          accuracy_order::Int, regfactor::Real) where {CM,K,DM}
        return new{CM,K,DM}(kernel, cm, dmgmodel, maxdmg, accuracy_order, regfactor)
    end
end

function RKCRMaterial(; kernel::Function=cubic_b_spline_kernel,
                        model::AbstractConstitutiveModel=LinearElastic(),
                        dmgmodel::AbstractDamageModel=CriticalStretch(), maxdmg::Real=0.85,
                        accuracy_order::Int=1, regfactor::Real=1e-13)
    if !(model isa LinearElastic)
        msg = "only `LinearElastic` is supported as model for `RKCRMaterial`!\n"
        throw(ArgumentError(msg))
    end
    if !(accuracy_order in (1, 2))
        msg = "RK kernel only implemented for `accuracy_order ∈ {1,2}`!\n"
        throw(ArgumentError(msg))
    end
    if !(0 ≤ regfactor ≤ 1)
        msg = "Regularization factor must be in the range 0 ≤ regfactor ≤ 1\n"
        throw(ArgumentError(msg))
    end
    return RKCRMaterial(kernel, model, dmgmodel, maxdmg, accuracy_order, regfactor)
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
    too_much_damage!(storage, system, mat, Fi, i) && return zero(SMatrix{3,3,Float64,9})
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
            Fb = 0.5 * (Fi + Fj)
            Ḟb = 0.5 * (Ḟi + Ḟj)
            # ΔXijLL = ΔXij' / (L * L)
            # ΔFij = (Δxij - Fb * ΔXij) * ΔXijLL
            # ΔḞij = (Δvij - Ḟb * ΔXij) * ΔXijLL
            # Fij = Fb + ΔFij
            # _Fij = Fb + ΔFij
            # Ḟij = Ḟb + ΔḞij
            # _Ḟij = Ḟb + ΔḞij

            LL = L * L

            β1 = Fb[1] * ΔXij[1] + Fb[2] * ΔXij[2] + Fb[3] * ΔXij[3]
            Fij1 = Fb[1] + (Δxij[1] - β1) * ΔXij[1] / LL
            Fij2 = Fb[2] + (Δxij[1] - β1) * ΔXij[2] / LL
            Fij3 = Fb[3] + (Δxij[1] - β1) * ΔXij[3] / LL

            β2 = Fb[4] * ΔXij[1] + Fb[5] * ΔXij[2] + Fb[6] * ΔXij[3]
            Fij4 = Fb[4] + (Δxij[2] - β2) * ΔXij[1] / LL
            Fij5 = Fb[5] + (Δxij[2] - β2) * ΔXij[2] / LL
            Fij6 = Fb[6] + (Δxij[2] - β2) * ΔXij[3] / LL

            β3 = Fb[7] * ΔXij[1] + Fb[8] * ΔXij[2] + Fb[9] * ΔXij[3]
            Fij7 = Fb[7] + (Δxij[3] - β3) * ΔXij[1] / LL
            Fij8 = Fb[8] + (Δxij[3] - β3) * ΔXij[2] / LL
            Fij9 = Fb[9] + (Δxij[3] - β3) * ΔXij[3] / LL

            Fij = SMatrix{3,3,Float64,9}(Fij1, Fij4, Fij7,
                                         Fij2, Fij5, Fij8,
                                         Fij3, Fij6, Fij9)

            β1 = Ḟb[1] * ΔXij[1] + Ḟb[2] * ΔXij[2] + Ḟb[3] * ΔXij[3]
            Ḟij1 = Ḟb[1] + (Δvij[1] - β1) * ΔXij[1] / LL
            Ḟij2 = Ḟb[2] + (Δvij[1] - β1) * ΔXij[2] / LL
            Ḟij3 = Ḟb[3] + (Δvij[1] - β1) * ΔXij[3] / LL

            β2 = Ḟb[4] * ΔXij[1] + Ḟb[5] * ΔXij[2] + Ḟb[6] * ΔXij[3]
            Ḟij4 = Ḟb[4] + (Δvij[2] - β2) * ΔXij[1] / LL
            Ḟij5 = Ḟb[5] + (Δvij[2] - β2) * ΔXij[2] / LL
            Ḟij6 = Ḟb[6] + (Δvij[2] - β2) * ΔXij[3] / LL

            β3 = Ḟb[7] * ΔXij[1] + Ḟb[8] * ΔXij[2] + Ḟb[9] * ΔXij[3]
            Ḟij7 = Ḟb[7] + (Δvij[3] - β3) * ΔXij[1] / LL
            Ḟij8 = Ḟb[8] + (Δvij[3] - β3) * ΔXij[2] / LL
            Ḟij9 = Ḟb[9] + (Δvij[3] - β3) * ΔXij[3] / LL

            Ḟij = SMatrix{3,3,Float64,9}(Ḟij1, Ḟij4, Ḟij7,
                                         Ḟij2, Ḟij5, Ḟij8,
                                         Ḟij3, Ḟij6, Ḟij9)

            # Fij ≈ _Fij || (@show Fij; @show _Fij)
            # Ḟij ≈ _Ḟij || (@show Ḟij; @show _Ḟij)

            # ORDERING OF THE ELEMENTS IN THE MATRICES IS IMPORTANT!
            #   // The default ordering is row-major:
            #   //
            #   // XX(0) XY(1) XZ(2)
            #   // YX(3) YY(4) YZ(5)
            #   // ZX(6) ZY(7) ZZ(8)

            # scalarTemp = *(meanDefGrad+0) * undeformedBondX + *(meanDefGrad+1) * undeformedBondY + *(meanDefGrad+2) * undeformedBondZ;
            # *(defGrad+0) = *(meanDefGrad+0) + (defStateX - scalarTemp) * undeformedBondX/undeformedBondLengthSq;
            # *(defGrad+1) = *(meanDefGrad+1) + (defStateX - scalarTemp) * undeformedBondY/undeformedBondLengthSq;
            # *(defGrad+2) = *(meanDefGrad+2) + (defStateX - scalarTemp) * undeformedBondZ/undeformedBondLengthSq;

            # scalarTemp = *(meanDefGrad+3) * undeformedBondX + *(meanDefGrad+4) * undeformedBondY + *(meanDefGrad+5) * undeformedBondZ;
            # *(defGrad+3) = *(meanDefGrad+3) + (defStateY - scalarTemp) * undeformedBondX/undeformedBondLengthSq;
            # *(defGrad+4) = *(meanDefGrad+4) + (defStateY - scalarTemp) * undeformedBondY/undeformedBondLengthSq;
            # *(defGrad+5) = *(meanDefGrad+5) + (defStateY - scalarTemp) * undeformedBondZ/undeformedBondLengthSq;

            # scalarTemp = *(meanDefGrad+6) * undeformedBondX + *(meanDefGrad+7) * undeformedBondY + *(meanDefGrad+8) * undeformedBondZ;
            # *(defGrad+6) = *(meanDefGrad+6) + (defStateZ - scalarTemp) * undeformedBondX/undeformedBondLengthSq;
            # *(defGrad+7) = *(meanDefGrad+7) + (defStateZ - scalarTemp) * undeformedBondY/undeformedBondLengthSq;
            # *(defGrad+8) = *(meanDefGrad+8) + (defStateZ - scalarTemp) * undeformedBondZ/undeformedBondLengthSq;

            # scalarTemp = *(meanDefGradDot+0) * undeformedBondX + *(meanDefGradDot+1) * undeformedBondY + *(meanDefGradDot+2) * undeformedBondZ;
            # *(defGradDot+0) = *(meanDefGradDot+0) + (velStateX - scalarTemp) * undeformedBondX/undeformedBondLengthSq;
            # *(defGradDot+1) = *(meanDefGradDot+1) + (velStateX - scalarTemp) * undeformedBondY/undeformedBondLengthSq;
            # *(defGradDot+2) = *(meanDefGradDot+2) + (velStateX - scalarTemp) * undeformedBondZ/undeformedBondLengthSq;

            # scalarTemp = *(meanDefGradDot+3) * undeformedBondX + *(meanDefGradDot+4) * undeformedBondY + *(meanDefGradDot+5) * undeformedBondZ;
            # *(defGradDot+3) = *(meanDefGradDot+3) + (velStateY - scalarTemp) * undeformedBondX/undeformedBondLengthSq;
            # *(defGradDot+4) = *(meanDefGradDot+4) + (velStateY - scalarTemp) * undeformedBondY/undeformedBondLengthSq;
            # *(defGradDot+5) = *(meanDefGradDot+5) + (velStateY - scalarTemp) * undeformedBondZ/undeformedBondLengthSq;

            # scalarTemp = *(meanDefGradDot+6) * undeformedBondX + *(meanDefGradDot+7) * undeformedBondY + *(meanDefGradDot+8) * undeformedBondZ;
            # *(defGradDot+6) = *(meanDefGradDot+6) + (velStateZ - scalarTemp) * undeformedBondX/undeformedBondLengthSq;
            # *(defGradDot+7) = *(meanDefGradDot+7) + (velStateZ - scalarTemp) * undeformedBondY/undeformedBondLengthSq;
            # *(defGradDot+8) = *(meanDefGradDot+8) + (velStateZ - scalarTemp) * undeformedBondZ/undeformedBondLengthSq;

            Pij = calc_first_piola_kirchhoff!(storage, mat, params, Fij, Ḟij, Δt, bond_id)
            # Tempij = I - ΔXij * ΔXijLL
            Tempij = I - (ΔXij * ΔXij') / LL
            wj = weighted_volume[j]
            ϕ = (wi > 0 && wj > 0) ? (0.5 / wi + 0.5 / wj) : 0.0
            # ϕ = (0.5 / wi + 0.5 / wj)
            ViVj = volume[i] * volume[j]
            ω̃ij = kernel(system, bond_id) * ϕ * ViVj
            ∑Pij = ω̃ij * (Pij * Tempij)
            ∑P += ∑Pij
            # @autoinfiltrate containsnan(∑P)
        end
    end
    return ∑P
end

function calc_first_piola_kirchhoff!(storage::RKCRStorage, mat::RKCRMaterial,
                                     params::StandardPointParameters, F::SMatrix{3,3,FT,9},
                                     Ḟ::SMatrix{3,3,FT,9}, Δt, bond_id) where {FT}
    D = init_stress_rotation!(storage, F, Ḟ, Δt, bond_id)
    iszero(D) && return zero(SMatrix{3,3,FT,9})
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
