# # [Writing your own material](@id tutorial_custom_material)

# This tutorial adds a material of your own to Peridynamics.jl: a bond-based material with
# a conical micro-modulus. Everything here is written against the [Extension API](@ref), so
# the same file runs single-threaded, with `julia -t 6`, and under
# `mpiexec -n 6 julia --project` without a change. The tutorials
# [Writing your own damage model](@ref tutorial_custom_damage_model) and
# [Writing your own constitutive model](@ref tutorial_custom_constitutive_model) build on
# the same ideas.
#
# The names of the extension API are not exported. Importing the ones a file uses keeps the
# code readable, and the methods the material defines are written as
# `Peridynamics.force_density_point!`, so that it is visible where the package is extended.

using Peridynamics
using Peridynamics: BondSystem, each_bond_idx, get_params, get_n_loc_points,
                    get_vector_diff, update_add_vector!, surface_correction_factor,
                    current_bond_length, bond_is_active, get_neighbor,
                    reference_bond_length, dims

# ## The material
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
# take the built-in [`CriticalStretch`](@ref), which breaks a bond when its stretch exceeds a
# critical value.

struct ConicalBBMaterial{Correction,DM} <: Peridynamics.AbstractBondSystemMaterial{Correction}
    dmgmodel::DM
end

function ConicalBBMaterial{C}(; dmgmodel=CriticalStretch()) where {C}
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
# `dmg_params`. Which parameters those are is decided by the damage model, not by the
# material. With `CriticalStretch` they are `Gc` and `εc`.

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

# Everything that can be inherited is listed in [Blocks you can inherit](@ref). For a
# material that is not restricted to ``\nu = 1/4``, `@inherit StandardParameters` gives the
# same set with the general elastic parameters, which take any two of `E`, `nu`, `G`, `K`,
# `λ` and `μ`, and it already includes the `dmg_params` marker.

# ### The storage
#
# The storage holds every field that changes during the simulation, one array per quantity.
# [`@storage`](@ref Peridynamics.@storage) generates the struct, its type parameters, the
# allocation, the halo exchange lists and the `Adapt` rule.
#
# `@inherit` pulls in the fields of the three time solvers, so the material works with all
# of them. The fracture bookkeeping is not declared here: it belongs to the damage model,
# and the `dmg_state::DamageState` marker is where the model puts it, whichever model the
# material is used with. We add one bond field of our own, so that the tutorial can also
# show how a field that is not a standard output reaches a VTK file.

Peridynamics.@storage ConicalBBMaterial struct ConicalBBStorage
    @inherit VelocityVerletFields DynamicRelaxationFields NewtonKrylovFields
    stretch::BondScalar
    dmg_state::DamageState
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
#
# What a force density may use is small and worth knowing by heart:
#
# - **The system**: `each_bond_idx(system, i)` iterates the bonds of a point,
#   `get_neighbor(system, bond_id)` is the point at the other end of a bond and
#   `reference_bond_length(system, bond_id)` its initial length. `system.volume[j]` is the
#   volume of a point, `system.position` the reference positions, `kernel(system, bond_id)`
#   the influence function and `surface_correction_factor(system, bond_id)` the surface
#   correction.
# - **The kinematics of a bond**: `current_bond_length(storage, system, i, bond_id)` is the
#   distance of the two points of the bond right now. Read it this way and never by gathering
#   the two positions yourself, then the same line is as fast as it can be on a material that
#   caches bond lengths and on one that does not. A kernel that needs the stretch and not the
#   length takes `bond_stretch(storage, system, i, bond_id)` instead, which is what the damage
#   model tutorial does.
# - **The storage**, through the fields of the blocks it inherited and its own fields, all
#   read the same flat way. Here these are `storage.position` and `storage.b_int` from the
#   solver blocks and `storage.stretch`. The table of every block says which fields it
#   brings. Whether a bond is broken is not read from a field but asked with
#   `bond_is_active(storage, system, bond_id)`, because the bookkeeping belongs to the
#   damage model and a material does not know which model it runs with.
# - **The parameters**: `paramsetup` is the parameter setup of the chunk, and
#   `get_params(paramsetup, i)` resolves the set of point `i`, with everything the `@params`
#   block declared.

function Peridynamics.force_density_point!(storage::ConicalBBStorage, system::BondSystem,
                                           mat::ConicalBBMaterial, paramsetup, t, Δt, i)
    (; volume) = system
    params = get_params(paramsetup, i)
    for bond_id in each_bond_idx(system, i)
        j = get_neighbor(system, bond_id)
        L = reference_bond_length(system, bond_id)

        ## the current bond vector, the current length and the stretch of the bond
        Δxij = get_vector_diff(storage.position, i, j, dims(system))
        l = current_bond_length(storage, system, i, bond_id)
        ε = (l - L) / L
        storage.stretch[bond_id] = ε

        ## the conical micro-modulus: full stiffness at L = 0, zero at the horizon
        c = params.bc * (1 - L / params.δ)

        ## a broken bond carries no force, and the surface correction is 1 for `NoCorrection`
        ω = bond_is_active(storage, system, bond_id) *
            surface_correction_factor(system, bond_id)

        ## the bond force, accumulated into point `i`
        b = ω * c * ε * volume[j] / l .* Δxij
        update_add_vector!(storage.b_int, i, b, dims(system))
    end
    return nothing
end

# Which bonds are broken was decided right before this call by the damage model, which is
# why the force density only multiplies `bond_is_active` in and never breaks anything.

# ### The fracture parameters
#
# A user gives a critical energy release rate ``G_c``, and `CriticalStretch` breaks a bond at
# a critical stretch ``\varepsilon_c``. The conversion between the two is an integral over
# all the bonds that a unit of crack surface cuts,
#
# ```math
# G_c = \frac{\pi \varepsilon_c^2}{2} \int_0^\delta c(\xi) \, \xi^4 \, \mathrm{d}\xi ,
# ```
#
# and it therefore depends on the micro-modulus. The package evaluates it for the constant
# micro-modulus by default and gets the familiar ``\varepsilon_c = \sqrt{5 G_c / (9 K \delta)}``.
# For the conical one, ``\int_0^\delta c_1 (1 - \xi/\delta) \xi^4 \mathrm{d}\xi = c_1 \delta^5/30``,
# so
#
# ```math
# G_c = \frac{3}{2} K \delta \varepsilon_c^2
# \qquad \Longleftrightarrow \qquad
# \varepsilon_c = \sqrt{\frac{2 G_c}{3 K \delta}} ,
# ```
#
# which is ``\sqrt{1.2} \approx 1.0954`` times the constant micro-modulus value. Put the
# other way round, the default conversion would return a critical stretch 8.7 % *below* the
# one this material implies, so the same ``G_c`` would break the body too early, silently.
#
# The relation follows from the micro-modulus, so the material states it, and it is the
# damage model that decides what ``\\varepsilon_c`` means, so both are dispatched on. Two
# methods, one for each direction, and `material!(...; Gc)` as well as
# `material!(...; epsilon_c)` convert correctly from now on. Define **both**: whichever of
# the two keywords the user gives, the other one is derived from it.

Peridynamics.critical_stretch(::CriticalStretch, ::ConicalBBMaterial, δ, K, Gc) =
    sqrt(2 * Gc / (3 * K * δ))

Peridynamics.energy_release_rate(::CriticalStretch, ::ConicalBBMaterial, δ, K, εc) =
    1.5 * K * δ * εc^2

# Nothing else about fracture has to be written. The damage model owns the keywords `Gc` and
# `epsilon_c`, decides which of them is an error to combine and switches fracture off when
# neither is given, all in the default
# [`get_frac_params`](@ref Peridynamics.get_frac_params). A method of that is only for a
# model that reads *other* keywords.

# ### Exporting a field of your own
#
# `stretch` is a bond field, so it cannot be written to a VTK file directly, because a
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
    weighted_stretch = zeros(get_n_loc_points(system))
    for i in eachindex(weighted_stretch)
        (; δ) = get_params(paramsetup, i)
        num, den = 0.0, 0.0
        for bond_id in each_bond_idx(system, i)
            L = reference_bond_length(system, bond_id)
            c = 1 - L / δ
            num += c * storage.stretch[bond_id]
            den += c
        end
        weighted_stretch[i] = iszero(den) ? 0.0 : num / den
    end
    return weighted_stretch
end

# ## Running it
#
# The material is finished. From here on nothing is specific to it. It is used exactly like
# a material that ships with the package.

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

# ## Where to go next
#
# - [Writing your own damage model](@ref tutorial_custom_damage_model) replaces the failure
#   criterion.
# - [Writing your own constitutive model](@ref tutorial_custom_constitutive_model) is the
#   way to go for a new stress-strain relation on the correspondence materials.
# - [Materials](@ref) for the declaration language in full.
# - [Blocks you can inherit](@ref) for everything `@inherit` accepts.
# - [Extension API](@ref) for every name used here.
# - [API stability](@ref) for what is promised and what is deliberately left internal.
