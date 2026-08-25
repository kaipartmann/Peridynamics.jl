# # [Writing your own material](@id tutorial_custom_material)

# This tutorial adds a material, a damage model and a constitutive model of your own to
# Peridynamics.jl. Everything here is written against the [Extension API](@ref), so the same
# file runs single-threaded, with `julia -t 6`, and under `mpiexec -n 6 julia --project`
# without a single change.

using Peridynamics
## `LinearAlgebra` and `StaticArrays` are reached through the package, so they do not have to
## be dependencies of your own project
using Peridynamics.LinearAlgebra

# ## Part 1: a custom material
#
# The bond-based formulation that [`BBMaterial`](@ref) implements uses a constant
# micro-modulus: every bond in the family of a point is equally stiff, no matter how long it
# is. That is the original choice of [Silling2005](@cite), and it is not the only one. A
# **conical** micro-modulus
#
# ```math
# c(\xi) = c_1 \left( 1 - \frac{\xi}{\delta} \right) , \qquad \xi = |\boldsymbol{\Delta X}| ,
# ```
#
# makes a bond stiff when it is short and lets its stiffness fall linearly to zero at the
# horizon, which removes the jump in stiffness at the edge of the family. It is one of the
# micro-modulus functions studied by [Bobaru2009](@cite) and [Ha2010](@cite), and it is what
# we build here.
#
# The constant ``c_1`` is not free. It follows from requiring that the material stores the
# same strain energy as a classical isotropic solid under a homogeneous stretch ``s``,
#
# ```math
# W = \pi s^2 \int_0^\delta c(\xi) \, \xi^3 \, \mathrm{d}\xi \overset{!}{=} \frac{9}{2} K s^2 .
# ```
#
# With ``\int_0^\delta c_1 (1 - \xi/\delta) \xi^3 \mathrm{d}\xi = c_1 \delta^4 / 20`` this
# gives ``c_1 = 90 K / (\pi \delta^4)``, five times the constant micro-modulus value
# ``18 K / (\pi \delta^4)``.
#
# A material needs four things: a type, its point parameters, its storage, and the force
# density calculation.

# ### The type
#
# The supertype decides which system the material is discretized on and therefore which
# arguments [`force_density_point!`](@ref Peridynamics.force_density_point!) gets. Bond-based
# means a bond system, so we subtype
# [`AbstractBondSystemMaterial`](@ref Peridynamics.AbstractBondSystemMaterial). Its type
# parameter is the surface correction, and we accept whatever the user asks for.
#
# The damage model is a field, because every bond system material is asked for one. We
# default it to the `ConicalStretch` model of Part 2, which is written for exactly this
# micro-modulus.

struct ConicalBBMaterial{Correction,DM} <: Peridynamics.AbstractBondSystemMaterial{Correction}
    dmgmodel::DM
end

function ConicalBBMaterial{C}(; dmgmodel=ConicalStretch()) where {C}
    return ConicalBBMaterial{C,typeof(dmgmodel)}(dmgmodel)
end
ConicalBBMaterial(; kwargs...) = ConicalBBMaterial{NoCorrection}(; kwargs...)

# ### The point parameters
#
# These are the values [`material!`](@ref) assigns to a point set.
# [`@params`](@ref Peridynamics.@params) generates the struct, the constructor that reads the
# keywords, the list of allowed keywords and the simulation log lines from one list of
# declarations, so they cannot drift apart.
#
# We compose two of the blocks the package ships: `DiscretizationParameters` gives `δ` and
# `rho`, and `BBElasticParameters` gives the six elastic constants in the ``\nu = 1/4``
# variant that bond-based peridynamics is restricted to. On top of that we derive the bond
# constant, and we give the parameters of the damage model a place with the marker field
# `dmg_params`. Which parameters those are is decided by the damage model in Part 2, not
# by the material.

Peridynamics.@params ConicalBBMaterial struct ConicalBBPointParameters
    @inherit DiscretizationParameters BBElasticParameters
    @log "conical micro-modulus constant" @derived bc = 90 * K / (π * δ^4)
    dmg_params::DamageParameters
end

# That is the whole parameter half. Three things happened in those five lines:
#
# - `material!(body; horizon, rho, E, Gc)` now works, and rejects anything misspelled. The
#   list of allowed keywords is derived from these declarations, so it cannot fall out of
#   sync with them.
# - `bc` is computed from parameters declared *before* it, and `@derived` keeps it out of
#   the keywords, because it is not something a user sets.
# - `@log` writes `bc` to the simulation log under the label given, so the value that was
#   actually used is recorded next to the results.
#
# The name `bc` is not free: every material on a bond system is asked for a bond constant
# under that name when the stable time step is estimated. Here `bc` is the micro-modulus at
# ``\xi = 0``, its largest value, so the estimate treats every bond as the stiffest one and
# returns a time step that is about ``\sqrt{3}`` smaller than it has to be. That is on the
# safe side, which is where a time step estimate belongs.
#
# Which names a right-hand side such as `90 * K / (π * δ^4)` may use is exactly what the
# inherited blocks expose, and every block says so itself:

Peridynamics.DiscretizationParameters

# The blocks this package ships are listed in [Blocks you can inherit](@ref). For a material
# that is not restricted to ``\nu = 1/4``, `@inherit StandardParameters` gives the same set
# with the general elastic parameters, which take any two of `E`, `nu`, `G`, `K`, `λ` and
# `μ`, and it already includes the `dmg_params` marker.

# ### The storage
#
# The storage holds every field that changes during the simulation, one array per quantity.
# [`@storage`](@ref Peridynamics.@storage) generates the struct, its type parameters, the
# allocation, the halo exchange lists and the `Adapt` rule.
#
# `@inherit` pulls in the fields of the three time solvers and of the fracture bookkeeping,
# so the material works with all of them. We add one bond field of our own, so that Part 1
# can also show how a field that is not a standard output reaches a VTK file.

Peridynamics.@storage ConicalBBMaterial struct ConicalBBStorage
    @inherit VelocityVerletFields DynamicRelaxationFields NewtonKrylovFields
    @inherit BondFracFields
    bond_stretch::BondScalar
end

# A field shape such as `BondScalar` says what a field means: how many entries it has, what
# its element type is, and how it is allocated. So no `init_field` method is needed. Without
# an explicit element type the field follows the float type of the simulation.

# ### The force density
#
# This is the one function a material has to define. It is called once per **local** point
# and per time step, inside a loop that runs on every thread and on every MPI rank. It only
# ever writes to the columns of point `i` and reads everything else through the system. That
# is the whole reason the same code parallelizes.

function Peridynamics.force_density_point!(storage::ConicalBBStorage,
                                           system::Peridynamics.BondSystem,
                                           mat::ConicalBBMaterial,
                                           params::ConicalBBPointParameters, t, Δt, i)
    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length

        ## the current bond vector and the bond stretch
        Δxij = Peridynamics.get_vector_diff(storage.position, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        storage.bond_stretch[bond_id] = ε

        ## the conical micro-modulus: full stiffness at L = 0, zero at the horizon
        c = params.bc * (1 - L / params.δ)

        ## a broken bond carries no force, and the surface correction is 1 for `NoCorrection`
        ω = storage.bond_active[bond_id] *
            Peridynamics.surface_correction_factor(system.correction, bond_id)

        ## the bond force, accumulated into point `i`
        b = ω * c * ε * system.volume[j] / l .* Δxij
        Peridynamics.update_add_vector!(storage.b_int, i, b)
    end
    return nothing
end

# ### Exporting a field of your own
#
# `bond_stretch` is a bond field, so it cannot be written to a VTK file directly, because a
# VTK file wants one value per point. [`export_field`](@ref Peridynamics.export_field)
# reduces it, and [`custom_field`](@ref Peridynamics.custom_field) announces the name so that
# asking for it in a `Job` is not rejected as a typo.
#
# We weight the reduction with the same micro-modulus the force uses, so that the exported
# number is the stretch the stiffness of the point actually sees, not a plain average over a
# family whose outer bonds barely carry load.

Peridynamics.custom_field(::Type{<:ConicalBBStorage}, ::Val{:weighted_stretch}) = true

function Peridynamics.export_field(::Val{:weighted_stretch}, mat, system,
                                   storage::ConicalBBStorage, paramsetup, t)
    weighted_stretch = zeros(Peridynamics.get_n_loc_points(system))
    for i in eachindex(weighted_stretch)
        δ = Peridynamics.get_params(paramsetup, i).δ
        num, den = 0.0, 0.0
        for bond_id in Peridynamics.each_bond_idx(system, i)
            c = 1 - system.bonds[bond_id].length / δ
            num += c * storage.bond_stretch[bond_id]
            den += c
        end
        weighted_stretch[i] = iszero(den) ? 0.0 : num / den
    end
    return weighted_stretch
end

# ## Part 2: a custom damage model
#
# A bond-based material breaks a bond when its stretch exceeds a critical value
# ``\varepsilon_c``, and the user would rather give a critical energy release rate ``G_c``
# than that stretch. The conversion between the two is an integral over all the bonds that
# a unit of crack surface cuts,
#
# ```math
# G_c = \frac{\pi \varepsilon_c^2}{2} \int_0^\delta c(\xi) \, \xi^4 \, \mathrm{d}\xi ,
# ```
#
# and it therefore depends on the micro-modulus. The built-in [`CriticalStretch`](@ref)
# evaluates it for the constant micro-modulus and gets the familiar
# ``\varepsilon_c = \sqrt{5 G_c / (9 K \delta)}``. For the conical one,
# ``\int_0^\delta c_1 (1 - \xi/\delta) \xi^4 \mathrm{d}\xi = c_1 \delta^5/30``, so
#
# ```math
# G_c = \frac{3}{2} K \delta \varepsilon_c^2
# \qquad \Longleftrightarrow \qquad
# \varepsilon_c = \sqrt{\frac{2 G_c}{3 K \delta}} ,
# ```
#
# which is ``\sqrt{1.2} \approx 1.0954`` times the constant micro-modulus value. Put the
# other way round, the built-in conversion returns a critical stretch 8.7 % *below* the one
# this material actually implies, so the same ``G_c`` would break the body too early. The
# material needs a damage model of its own.
#
# A damage model is a type, its parameters and the conversion. The parameters are the
# standard fracture parameters `Gc` and `εc`, which the
# [`FractureParameters`](@ref Peridynamics.FractureParameters) block declares. It reads the
# keywords `Gc` and `epsilon_c` of `material!` and hands them to
# [`get_frac_params`](@ref Peridynamics.get_frac_params) of the model, which is where the
# conversion lives. The failure criterion itself is unchanged, a bond breaks when it is
# stretched too far, so `calc_failure!` and `has_fracture` forward to the built-in model
# rather than repeating it.

struct ConicalStretch <: Peridynamics.AbstractDamageModel end

Peridynamics.@dmg_params ConicalStretch struct ConicalStretchParameters
    @inherit FractureParameters
end

function Peridynamics.get_frac_params(::ConicalStretch, δ, K; Gc=nothing, epsilon_c=nothing,
                                      kwargs...)
    if !isnothing(Gc) && isnothing(epsilon_c)
        return (; Gc=float(Gc), εc=sqrt(2 * Gc / (3 * K * δ)))
    elseif isnothing(Gc) && !isnothing(epsilon_c)
        return (; Gc=1.5 * K * δ * epsilon_c^2, εc=float(epsilon_c))
    elseif !isnothing(Gc) && !isnothing(epsilon_c)
        throw(ArgumentError("define either Gc or epsilon_c, not both!\n"))
    end
    return (; Gc=0.0, εc=0.0)
end

function Peridynamics.calc_failure!(storage, system, mat, ::ConicalStretch, paramsetup, i)
    return Peridynamics.calc_failure!(storage, system, mat, CriticalStretch(), paramsetup, i)
end

function Peridynamics.has_fracture(::ConicalStretch, params)
    return Peridynamics.has_fracture(CriticalStretch(), params)
end

# `material!` now accepts `Gc` or `epsilon_c` for a `ConicalBBMaterial` and rejects both at
# once. The keywords a damage model reads are the keyword arguments of its own
# `get_frac_params` method, and a keyword the user did not give arrives as `nothing`.

# ### Running it
#
# The material and its damage model are finished. From here on nothing is specific to them.
# They are used exactly like a material that ships with the package.

l, Δx = 0.1, 0.002
pos, vol = uniform_box(l, 0.1l, 0.1l, Δx)
body = Body(ConicalBBMaterial(), pos, vol)
material!(body; horizon=3.015Δx, rho=2700, E=70e9, Gc=100)

# Pull the bar apart at both ends:

point_set!(x -> x < -0.4l, body, :left)
point_set!(x -> x > 0.4l, body, :right)
velocity_bc!(t -> -10.0, body, :left, :x)
velocity_bc!(t -> 10.0, body, :right, :x)

# Ask for our own field alongside the built-in ones:

job = Job(body, VelocityVerlet(steps=200);
          path=joinpath(tempdir(), "custom_material"),
          fields=(:displacement, :damage, :weighted_stretch))

#-
#md # ```julia
#md # submit(job)
#md # ```

# The same file, unchanged, also runs as
#
# ```bash
# julia -t 6 --project tutorial_custom_material.jl        # multithreading
# mpiexec -n 6 julia --project tutorial_custom_material.jl # MPI
# ```
#
# You did not write anything for that. The halo exchange follows from the field
# declarations: `position` is annotated `@lth` in `VelocityVerletFields`, so it is sent from
# the chunk that owns a point to the chunks that need it, and everything else is chunk-local.

# ## Part 3: a custom constitutive model
#
# For the correspondence families ([`CMaterial`](@ref), [`RKCMaterial`](@ref),
# [`BACMaterial`](@ref)) a new stress-strain relation usually does **not** need a new
# material at all. Those materials ask a constitutive model for the first Piola-Kirchhoff
# stress that belongs to a deformation gradient, so the model is all you write, and it then
# runs on every one of those families.
#
# As an example, the Mooney-Rivlin solid [Mooney1940](@cite), [Rivlin1948](@cite), in the
# volumetric-isochoric split
#
# ```math
# W = C_{10} (\bar{I}_1 - 3) + C_{01} (\bar{I}_2 - 3) + \frac{K}{2} (J - 1)^2 ,
# ```
#
# whose second Piola-Kirchhoff stress is ``\boldsymbol{S} = 2 \, \partial W / \partial
# \boldsymbol{C}``:

struct MooneyRivlin <: Peridynamics.AbstractConstitutiveModel
    C10::Float64
    C01::Float64
end
MooneyRivlin(; C10=0.3, C01=0.1) = MooneyRivlin(C10, C01)

function Peridynamics.first_piola_kirchhoff(cm::MooneyRivlin, storage, params, F)
    J = det(F)
    J < eps() && return zero(F)
    Finv = inv(F)
    ## the first two invariants of the right Cauchy-Green tensor
    C = F' * F
    I1, I2 = tr(C), 0.5 * (tr(C)^2 - tr(C * C))
    P = 2 * cm.C10 * J^(-2 / 3) * (F - I1 / 3 * Finv') +
        2 * cm.C01 * J^(-4 / 3) * (I1 * F - F * C - 2I2 / 3 * Finv') +
        params.K * (J - 1) * J * Finv'
    return P
end

# `RKCMaterial(model=MooneyRivlin())` now works, and so does `CMaterial(model=…)` and
# `BACMaterial(model=…)`:

hyperelastic_body = Body(RKCMaterial(; model=MooneyRivlin()), pos, vol)
material!(hyperelastic_body; horizon=3.015Δx, rho=2700, E=70e9, nu=0.3, Gc=100)

# ### A model with a history
#
# A model that integrates an internal state over time, such as plasticity, viscoelasticity
# or creep, declares that state with [`@cm_storage`](@ref Peridynamics.@cm_storage), which
# takes the same field declarations as `@storage`:
#
# ```julia
# struct J2Plasticity <: Peridynamics.AbstractConstitutiveModel end
#
# Peridynamics.@cm_params J2Plasticity struct J2PlasticityParameters
#     @log "yield stress" sigma_y
#     @log "hardening modulus" H = 0.0
# end
#
# Peridynamics.@cm_storage J2Plasticity struct J2PlasticityState
#     bond_plastic_strain::BondSymTensor
#     bond_eqps::BondScalar
# end
# ```
#
# The parameters are declared with [`@cm_params`](@ref Peridynamics.@cm_params) and become
# keywords of `material!`. The state is allocated with the storage, moves with it to another
# array backend, and is reached inside the stress update, which then takes two more
# arguments: the index of the evaluated quantity and the time step. A symmetric field has
# six rows rather than nine, so it is read and written with
# [`get_sym_tensor`](@ref Peridynamics.get_sym_tensor) and
# [`update_sym_tensor!`](@ref Peridynamics.update_sym_tensor!).
#
# ```julia
# function Peridynamics.first_piola_kirchhoff(cm::J2Plasticity, storage, params, F, idx, Δt)
#     state = Peridynamics.constitutive_state(storage)
#     εᵖ = Peridynamics.get_sym_tensor(state.bond_plastic_strain, idx)
#     ## the elastic predictor in logarithmic strain space
#     ε, Uinv = Peridynamics.hencky_and_invstretch(F' * F)
#     τ = params.λ * tr(ε - εᵖ) * I + 2 * params.μ * (ε - εᵖ)
#     ## ... radial return with `params.sigma_y` and `params.H`, which updates `εᵖ` ...
#     Peridynamics.update_sym_tensor!(state.bond_plastic_strain, idx, εᵖ_new)
#     return F * (Uinv * τ * Uinv)
# end
# ```
#
# Declaring a state is what makes the model history dependent, and that is checked when the
# [`Job`](@ref) is created: a solver that evaluates the force density more than once per
# step, such as [`NewtonKrylov`](@ref), is rejected rather than integrating the history
# several times.

# ## Where to go next
#
# - [Materials](@ref) for the declaration language in full.
# - [Blocks you can inherit](@ref) for everything `@inherit` accepts.
# - [Extension API](@ref) for every name used here.
# - [API stability](@ref) for what is promised and what is deliberately left internal.
