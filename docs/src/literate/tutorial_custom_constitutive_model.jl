# # [Writing your own constitutive model](@id tutorial_custom_constitutive_model)

# For the correspondence materials ([`CMaterial`](@ref), [`RKCMaterial`](@ref),
# [`BACMaterial`](@ref)) a new stress-strain relation does not need a new material at all.
# Those materials ask a **constitutive model** for the first Piola-Kirchhoff stress that
# belongs to a deformation gradient, so the model is all you write, and it then runs on every
# one of those families, on threads and with MPI. This tutorial writes two: a hyperelastic
# model without state, and a plasticity model that carries a history.

using Peridynamics
using Peridynamics: constitutive_state, get_sym_tensor, update_sym_tensor!,
                    hencky_and_invstretch, each_bond_idx, get_n_loc_points
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
# has two constants of its own, ``C_{10}`` and ``C_{01}``, and reads the bulk modulus ``K``
# of the material. A constitutive model is a type, so the constants could be fields of the
# struct. Declaring them with [`@cm_params`](@ref Peridynamics.@cm_params) is better: they
# become keywords of [`material!`](@ref), they are written to the simulation log, and they
# can differ between point sets like every other material parameter.

struct MooneyRivlin <: Peridynamics.AbstractConstitutiveModel end

Peridynamics.@cm_params MooneyRivlin struct MooneyRivlinParameters
    @log "Mooney-Rivlin constant C10" C10
    @log "Mooney-Rivlin constant C01" C01 = 0.0
end

# `C10` is required and `C01` defaults to zero, which makes the model a neo-Hookean solid.
# The parameters are read flat off the point parameters, `params.C10`, next to the elastic
# constants the material declares, because the correspondence materials carry the marker
# field `cm_params::ConstitutiveParameters` for exactly this purpose.
#
# The stress follows from ``\boldsymbol{S} = 2 \, \partial W / \partial \boldsymbol{C}`` and
# ``\boldsymbol{P} = \boldsymbol{F} \boldsymbol{S}``. The one method a constitutive model
# has to define is [`first_piola_kirchhoff`](@ref Peridynamics.first_piola_kirchhoff), and
# a model without state defines the four-argument form:

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
# A model that integrates an internal state over time, such as plasticity, viscoelasticity
# or creep, needs three things more: parameters of its own, a state per evaluation point,
# and the time step. We write rate-independent J2 plasticity with linear isotropic
# hardening, integrated in logarithmic strain space with an elastic predictor and a radial
# return.

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
# The plastic strain and the equivalent plastic strain have to be remembered between time
# steps. [`@cm_storage`](@ref Peridynamics.@cm_storage) declares them with the field shapes
# of [`@storage`](@ref Peridynamics.@storage). The state is allocated with the storage of the
# material, moves with it to another array backend, and is never exchanged between chunks,
# because plastic state is local to the point it belongs to.
#
# The reproducing kernel materials evaluate the constitutive model once per bond, at the
# bond-associated quadrature point, so the state is declared per bond. A symmetric tensor
# has six independent components, and `BondSymTensor` stores exactly those.

Peridynamics.@cm_storage J2Plasticity struct J2PlasticityState
    bond_plastic_strain::BondSymTensor
    bond_eqps::BondScalar
end

# Declaring a state is what makes the model history dependent. That is checked when the
# [`Job`](@ref) is created: a time solver that evaluates the force density more than once
# per step, such as [`NewtonKrylov`](@ref), is rejected with a
# [`HistoryDependenceError`](@ref Peridynamics.HistoryDependenceError) rather than
# integrating the history several times.

# ### The stress update
#
# A model with a state defines the six-argument form of `first_piola_kirchhoff`, which also
# receives the index of the evaluated quantity and the time step. The state is reached with
# [`constitutive_state`](@ref Peridynamics.constitutive_state), and a symmetric field is read
# and written with [`get_sym_tensor`](@ref Peridynamics.get_sym_tensor) and
# [`update_sym_tensor!`](@ref Peridynamics.update_sym_tensor!).
#
# The elastic predictor works in logarithmic strain space:
# [`hencky_and_invstretch`](@ref Peridynamics.hencky_and_invstretch) returns the Hencky
# strain ``\boldsymbol{\varepsilon} = \ln \boldsymbol{U}`` and the inverse stretch
# ``\boldsymbol{U}^{-1}`` in closed form. The trial stress is the rotated Kirchhoff stress
# that belongs to the elastic part of the strain, so the yield surface caps the true stress.
# If the trial stress lies outside the yield surface, the radial return brings it back, and
# for linear hardening the plastic multiplier is known in closed form. The pull back to the
# first Piola-Kirchhoff stress is
# ``\boldsymbol{P} = \boldsymbol{F} \boldsymbol{U}^{-1} \hat{\boldsymbol{\tau}} \boldsymbol{U}^{-1}``.

function Peridynamics.first_piola_kirchhoff(::J2Plasticity, storage, params, F, idx, Δt)
    (; λ, μ, sigma_y, H) = params
    state = constitutive_state(storage)

    ## elastic predictor in logarithmic strain space
    ε, Uinv = hencky_and_invstretch(F' * F)
    εp = get_sym_tensor(state.bond_plastic_strain, idx)
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
    update_sym_tensor!(state.bond_plastic_strain, idx, εp + Δεp)
    state.bond_eqps[idx] = eqps + Δγ
    return F * (Uinv * τ * Uinv)
end

# The elastic parameters `λ` and `μ` come from the material, `sigma_y` and `H` from the
# model, and all four are read the same way. The strain energy density is what the export
# of the `:strain_energy_density` field asks for. It takes the index but not the time step,
# and it must not change the state:

function Peridynamics.strain_energy_density(::J2Plasticity, storage, params, F, idx)
    state = constitutive_state(storage)
    ε, _ = hencky_and_invstretch(F' * F)
    εe = ε - get_sym_tensor(state.bond_plastic_strain, idx)
    return 0.5 * params.λ * tr(εe)^2 + params.μ * tr(εe * εe)
end

# ### Exporting the plastic strain
#
# The equivalent plastic strain lives per bond. To see it in ParaView it is reduced to one
# value per point, here the maximum over the bonds of a point, because the worst bond is what
# governs whether a point is about to fail.

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

# Nothing in the model mentions a thread, a rank or the material family that evaluates it.
# Since the state is declared per bond, the model runs on the materials that evaluate their
# constitutive model per bond, [`RKCMaterial`](@ref) and [`BACMaterial`](@ref). For
# [`CMaterial`](@ref), which evaluates it once per point, the state would be declared with
# the point shapes instead.

# ## Where to go next
#
# - [Writing your own material](@ref tutorial_custom_material) for a material that is not a
#   constitutive model.
# - [Writing your own damage model](@ref tutorial_custom_damage_model) for a failure
#   criterion of your own, which can read the state of a constitutive model.
# - [Materials](@ref) for the declaration language in full.
# - [Extension API](@ref) for every name used here.
