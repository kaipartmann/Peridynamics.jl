# # [Writing your own constitutive model](@id tutorial_custom_constitutive_model)

# For the correspondence materials ([`CMaterial`](@ref), [`RKCMaterial`](@ref),
# [`BACMaterial`](@ref)) a new stress-strain relation needs no new material. These ask a
# **constitutive model** for the first Piola-Kirchhoff stress of a deformation gradient, and
# it runs on threads and MPI alike. This tutorial writes two: one without state, one with.

using Peridynamics
using Peridynamics: constitutive_state, get_sym_tensor, update_sym_tensor!,
                    hencky_and_invstretch, each_bond_idx, get_n_loc_points, dims
using Peridynamics.LinearAlgebra: norm, tr, det, inv, I

# ## A hyperelastic model
#
# The Mooney-Rivlin solid [Mooney1940](@cite), [Rivlin1948](@cite) in the
# volumetric-isochoric split
#
# ```math
# W = C_{10} (\bar{I}_1 - 3) + C_{01} (\bar{I}_2 - 3) + \frac{K}{2} (J - 1)^2
# ```
#
# has two constants of its own, ``C_{10}`` and ``C_{01}``, and reads the bulk modulus ``K``.
# Declaring them with [`@cm_params`](@ref Peridynamics.@cm_params) makes them keywords of
# [`material!`](@ref), logged, and free to differ between point sets.

struct MooneyRivlin <: Peridynamics.AbstractConstitutiveModel end

Peridynamics.@cm_params MooneyRivlin struct MooneyRivlinParameters
    @log "Mooney-Rivlin constant C10" C10
    @log "Mooney-Rivlin constant C01" C01 = 0.0
end

# `C10` is required and `C01` defaults to zero, a neo-Hookean solid. Parameters are read flat
# as `params.C10`, via the marker field `cm_params::ConstitutiveParameters`.
#
# The stress follows from ``\boldsymbol{S} = 2 \, \partial W / \partial \boldsymbol{C}`` and
# ``\boldsymbol{P} = \boldsymbol{F} \boldsymbol{S}``. The one method to define is
# [`first_piola_kirchhoff`](@ref Peridynamics.first_piola_kirchhoff), here its four-argument
# form for a model without state:

function Peridynamics.first_piola_kirchhoff(::MooneyRivlin, storage, params, F)
    J = det(F)
    J < eps() && return zero(F)
    Finv = inv(F)
    ## the first two invariants of the right Cauchy-Green tensor
    C = F' * F
    I1, I2 = tr(C), 0.5 * (tr(C)^2 - tr(C * C))
    P = 2 * params.C10 * J^(-2 / 3) * (F - I1 / 3 * Finv') +
        2 * params.C01 * J^(-4 / 3) * (I1 * F - F * C - 2I2 / 3 * Finv') +
        params.K * (J - 1) * J * Finv'
    return P
end

# That is the whole model. `RKCMaterial(model=MooneyRivlin())` works, and so does
# `CMaterial(model=MooneyRivlin())` and `BACMaterial(model=MooneyRivlin())`:

l, Δx = 0.1, 0.004
pos, vol = uniform_box(l, 0.2l, 0.2l, Δx)
rubber = Body(RKCMaterial(; model=MooneyRivlin()), pos, vol)
material!(rubber; horizon=3.015Δx, rho=1100, E=2.4e6, nu=0.49, C10=0.3e6, C01=0.1e6)

# ## A model with a history
#
# A model that integrates an internal state over time, such as plasticity or viscoelasticity,
# needs parameters of its own, a state per evaluation point, and the time step. We write
# rate-independent J2 plasticity with linear isotropic hardening and a radial return.

struct J2Plasticity <: Peridynamics.AbstractConstitutiveModel end

# ### The parameters
#
# The yield stress is required, the hardening modulus defaults to zero, which is perfect
# plasticity. Both become keywords of `material!`.

Peridynamics.@cm_params J2Plasticity struct J2PlasticityParameters
    @log "yield stress" sigma_y
    @log "hardening modulus" H = 0.0
end

# ### The state
#
# The plastic strain and equivalent plastic strain persist between time steps, declared with
# [`@cm_storage`](@ref Peridynamics.@cm_storage) using the field shapes of
# [`@storage`](@ref Peridynamics.@storage), and are never exchanged between chunks.
#
# The reproducing kernel materials evaluate the model once per bond, so the state is declared
# per bond, with `BondSymTensor` storing its six independent components.

Peridynamics.@cm_storage J2Plasticity struct J2PlasticityState
    bond_plastic_strain::BondSymTensor
    bond_eqps::BondScalar
end

# Declaring a state makes the model history dependent: a time solver evaluating the force
# density more than once per step, such as [`NewtonKrylov`](@ref), is rejected with a
# [`HistoryDependenceError`](@ref Peridynamics.HistoryDependenceError).

# ### The stress update
#
# A model with a state defines the six-argument form of `first_piola_kirchhoff`, receiving
# also the index and the time step. The state is reached with
# [`constitutive_state`](@ref Peridynamics.constitutive_state), read and written with
# [`get_sym_tensor`](@ref Peridynamics.get_sym_tensor) and
# [`update_sym_tensor!`](@ref Peridynamics.update_sym_tensor!).
#
# [`hencky_and_invstretch`](@ref Peridynamics.hencky_and_invstretch) gives the elastic
# predictor in logarithmic strain space. The trial stress, capped by the yield surface, gets
# a closed-form radial return for linear hardening when it lies outside it.

function Peridynamics.first_piola_kirchhoff(::J2Plasticity, storage, params, F, idx, Δt)
    (; λ, μ, sigma_y, H) = params
    state = constitutive_state(storage)

    ## elastic predictor in logarithmic strain space
    ε, Uinv = hencky_and_invstretch(F' * F)
    ## the constitutive model API has no system in scope, and this tutorial is 3D
    εp = get_sym_tensor(state.bond_plastic_strain, idx, dims(storage))
    εe = ε - εp
    τ_trial = λ * tr(εe) * I + 2 * μ * εe

    ## the deviatoric part decides about yielding, the pressure is unaffected
    p = tr(τ_trial) / 3
    s = τ_trial - p * I
    q = sqrt(1.5) * norm(s)
    eqps = state.bond_eqps[idx]
    σy = sigma_y + H * eqps
    q ≤ σy && return F * (Uinv * τ_trial * Uinv)

    ## radial return: the plastic multiplier is closed form for linear hardening
    Δγ = (q - σy) / (3 * μ + H)
    Δεp = sqrt(1.5) * Δγ * (s / norm(s))
    τ = τ_trial - 2 * μ * Δεp
    update_sym_tensor!(state.bond_plastic_strain, idx, εp + Δεp, dims(storage))
    state.bond_eqps[idx] = eqps + Δγ
    return F * (Uinv * τ * Uinv)
end

# `λ`, `μ`, `sigma_y` and `H` all read the same way. `:strain_energy_density` export takes
# the index, not the time step, and must not change the state:

function Peridynamics.strain_energy_density(::J2Plasticity, storage, params, F, idx)
    state = constitutive_state(storage)
    ε, _ = hencky_and_invstretch(F' * F)
    εe = ε - get_sym_tensor(state.bond_plastic_strain, idx, dims(storage))
    return 0.5 * params.λ * tr(εe)^2 + params.μ * tr(εe * εe)
end

# ### Exporting the plastic strain
#
# The equivalent plastic strain lives per bond, reduced to one value per point, the maximum
# over its bonds, since the worst bond governs whether a point is about to fail.

Peridynamics.custom_field(::Type{<:Peridynamics.AbstractStorage},
                          ::Val{:equivalent_plastic_strain}) = true

function Peridynamics.export_field(::Val{:equivalent_plastic_strain}, mat, system, storage,
                                   paramsetup, t)
    eqps = constitutive_state(storage).bond_eqps
    out = zeros(get_n_loc_points(system))
    for i in eachindex(out)
        bond_ids = each_bond_idx(system, i)
        isempty(bond_ids) || (out[i] = maximum(@view eqps[bond_ids]))
    end
    return out
end

# ## Running it
#
# A bar of steel is pulled far beyond its yield stress:

steel = Body(RKCMaterial(; model=J2Plasticity()), pos, vol)
material!(steel; horizon=3.015Δx, rho=7850, E=210e9, nu=0.3, sigma_y=250e6, H=1e9)
point_set!(x -> x < -0.4l, steel, :left)
point_set!(x -> x > 0.4l, steel, :right)
velocity_bc!(t -> -5.0, steel, :left, :x)
velocity_bc!(t -> 5.0, steel, :right, :x)
job = Job(steel, VelocityVerlet(steps=2000);
          path=joinpath(tempdir(), "plasticity"),
          fields=(:displacement, :damage, :equivalent_plastic_strain))

#-
#md # ```julia
#md # submit(job)
#md # ```

# Nothing in the model mentions a thread, a rank or the material family evaluating it. The
# state declared per bond runs on [`RKCMaterial`](@ref) and [`BACMaterial`](@ref).
# [`CMaterial`](@ref) would need it declared with the point shapes instead.

# ## Where to go next
#
# - [Writing your own material](@ref tutorial_custom_material) for a material that is not a
#   constitutive model.
# - [Writing your own damage model](@ref tutorial_custom_damage_model) for a failure
#   criterion of your own, which can read the state of a constitutive model.
# - [Materials](@ref) for the declaration language in full.
# - [Extension API](@ref) for every name used here.
