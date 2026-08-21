#=
The constitutive model interface.

A constitutive model answers one question: which stress belongs to a deformation gradient.
Everything a *history-dependent* model needs on top of that lives here as well:

- the canonical signatures, which carry the index of the evaluation point and the time step,
  so that a model can read and advance its own state,
- `@cm_storage` (`core/storages.jl`), with which a model declares that state the same way a
  material declares its storage,
- the traits `is_history_dependent` and `supports_history_dependence`, which are checked once
  when a `Job` is created, so that a model that integrates its history cannot silently be
  run by a solver that evaluates the force density several times per step.

The concrete models of this package are in `physics/constitutive_models.jl`. This file is
included early, so that a material family can declare the state of its constitutive model
next to the family itself.
=#

"""
    get_constitutive_model(mat)

$(extension_api_note())

Return the constitutive model of a material, or `nothing` for a material family that does
not have one, e.g. the bond-based materials. Every family that stores a constitutive model
defines this method, which is how the setup checks and the storage of the constitutive state
reach the model without knowing the material.

# Example
```julia
Peridynamics.get_constitutive_model(mat::MyMaterial) = mat.constitutive_model
```
"""
function get_constitutive_model end

get_constitutive_model(::AbstractMaterial) = nothing

# the marker has to be interpolated, the body has to be raw for the LaTeX, so the two halves
# are concatenated
@doc """
    first_piola_kirchhoff(model, storage, params, F, idx, Δt)

$(extension_api_note())
""" * raw"""
Return the first Piola-Kirchhoff stress ``\boldsymbol{P}`` of a constitutive model for the
deformation gradient ``\boldsymbol{F}``. This is the one method a constitutive model has to
define.

# Arguments

- `model`: The constitutive model.
- `storage`: The storage of the body chunk. The state of a history-dependent model is
    reached with [`constitutive_state`](@ref).
- `params`: The point parameters of the point that is evaluated.
- `F::SMatrix{3,3}`: The deformation gradient.
- `idx::Int`: The index of the evaluated quantity, which is what a history-dependent model
    indexes its state with. Which index this is follows from the material family:

    | material family | `idx` |
    |:---|:---|
    | [`CMaterial`](@ref) | point index |
    | [`RKCMaterial`](@ref), [`BACMaterial`](@ref) | bond index |

- `Δt::Real`: The time step, i.e. the increment a history-dependent model integrates over.

# A model without state

A model that only needs the deformation gradient defines the four-argument form and ignores
this one, which is bridged to it:

```julia
function Peridynamics.first_piola_kirchhoff(::MyModel, storage, params, F)
    return ...
end
```

# A model with state

```julia
struct J2Plasticity <: Peridynamics.AbstractConstitutiveModel end

Peridynamics.@cm_storage J2Plasticity struct J2PlasticityState
    bond_plastic_strain::BondTensor
    bond_eqps::BondScalar
end

function Peridynamics.first_piola_kirchhoff(cm::J2Plasticity, storage, params, F, idx, Δt)
    state = Peridynamics.constitutive_state(storage)
    εᵖ = Peridynamics.get_tensor(state.bond_plastic_strain, idx)
    # ... radial return ...
    Peridynamics.update_tensor!(state.bond_plastic_strain, idx, εᵖ_new)
    return P
end
```
"""
function first_piola_kirchhoff end

# The bridge that keeps every stateless model, in this package and out of tree, working with
# the four-argument form it was written with. It is `@inline`, so the two extra arguments
# disappear completely.
@inline function first_piola_kirchhoff(model::AbstractConstitutiveModel, storage, params, F,
                                       idx, Δt)
    return first_piola_kirchhoff(model, storage, params, F)
end

function first_piola_kirchhoff(model::AbstractConstitutiveModel, storage, params, F)
    return throw(InterfaceError(model, "first_piola_kirchhoff",
                                constitutive_model_hint(model, "first_piola_kirchhoff")))
end

@doc """
    strain_energy_density(model, storage, params, F, idx)

$(extension_api_note())
""" * raw"""
Return the strain energy density ``\Psi`` of a constitutive model. In contrast to
[`first_piola_kirchhoff`](@ref) this method must not change the state of the model: it is
also called when a field is exported, i.e. outside of the time integration. See
[`first_piola_kirchhoff`](@ref) for the meaning of `idx`.

A model that only needs the deformation gradient defines the four-argument form
`strain_energy_density(model, storage, params, F)`, which this one is bridged to.
"""
function strain_energy_density end

@inline function strain_energy_density(model::AbstractConstitutiveModel, storage, params, F,
                                       idx)
    return strain_energy_density(model, storage, params, F)
end

function strain_energy_density(model::AbstractConstitutiveModel, storage, params, F)
    return throw(InterfaceError(model, "strain_energy_density",
                                constitutive_model_hint(model, "strain_energy_density")))
end

function constitutive_model_hint(model, func::AbstractString)
    name = model isa Type ? nameof(model) : nameof(typeof(model))
    msg = "Define the constitutive model `$(name)`, e.g.\n"
    msg *= "        function Peridynamics.$(func)(::$(name), storage, params, F)\n"
    msg *= "            return ...\n"
    msg *= "        end\n"
    msg *= "    or, for a history-dependent model, the form that also takes the index of "
    msg *= "the\n    evaluated quantity and the time step, see `?Peridynamics."
    msg *= "first_piola_kirchhoff`."
    return msg
end

# --------------------------------------------------------------------------------------
# state of a history-dependent model
# --------------------------------------------------------------------------------------

"""
    constitutive_storage_type(model)
    constitutive_storage_type(model, ::Type{FT})

$(extension_api_note())

Return the type of the state of a constitutive model, instantiated for the float type `FT`
of the simulation, or `Nothing` for a model without state. This method is generated by
[`@cm_storage`](@ref) and is the constitutive-model analogue of [`storage_type`](@ref).
"""
function constitutive_storage_type end

constitutive_storage_type(model) = Nothing
constitutive_storage_type(model, ::Type) = constitutive_storage_type(model)

"""
    get_cm_storage(model, solver, system)

$(internal_api_warning())

Allocate the state of a constitutive model for one body chunk, or return `nothing` for a
model without state. This method is generated by [`@cm_storage`](@ref) and is the
constitutive-model analogue of `get_storage`.
"""
function get_cm_storage end

get_cm_storage(model, solver, system) = nothing

"""
    init_constitutive_state(mat, solver, system)

$(internal_api_warning())

Return the value of the storage field that was declared with
`cm_state::ConstitutiveState`, which by default is the state of the constitutive model of
`mat`. This is the escape hatch for a material family that assembles the state of its
constitutive model itself.
"""
function init_constitutive_state(mat, solver, system)
    return get_cm_storage(get_constitutive_model(mat), solver, system)
end

"""
    constitutive_state(storage)

$(extension_api_note())

Return the state of the constitutive model that is carried by a storage, i.e. the field
declared with `cm_state::ConstitutiveState`, or `nothing` for a storage that does not
declare one. This is how [`first_piola_kirchhoff`](@ref) reaches the state of a
history-dependent model.
"""
function constitutive_state end

@inline constitutive_state(::AbstractStorage) = nothing

"""
    has_constitutive_state(::Type{Storage})

$(internal_api_warning())

Return whether a storage type declares a field with the [`ConstitutiveState`](@ref)
declaration and can therefore carry the state of a history-dependent constitutive model.
"""
function has_constitutive_state end

has_constitutive_state(::Type) = false

"""
    is_history_dependent(model)

$(extension_api_note())

Return whether a constitutive model integrates an internal state over time, e.g. the plastic
strain of a plasticity model. Defaults to whether the model declares a state with
[`@cm_storage`](@ref), so a model that keeps its history in the storage of the material has
to declare it:

```julia
Peridynamics.is_history_dependent(::MyPlasticModel) = true
```

A history-dependent model may only be used with a time solver that evaluates the force
density once per time step, and only with a material whose storage carries the state. Both
are checked once when a [`Job`](@ref) is created.
"""
function is_history_dependent end

is_history_dependent(model) = constitutive_storage_type(model) !== Nothing
is_history_dependent(::Nothing) = false

"""
    supports_history_dependence(solver)
    supports_history_dependence(mat)

$(extension_api_note())

Return whether a time solver or a material may be combined with a history-dependent
constitutive model. Defaults to `true` for both.

A time solver answers `false` if it evaluates the force density more than once per time
step, because the model would then integrate its history several times per step.
[`NewtonKrylov`](@ref) does this for the Jacobian-vector products and the line search.

A material answers `false` if the index that [`first_piola_kirchhoff`](@ref) receives does
not always index the same quantity, because the state of the model could then not be
indexed at all.
"""
function supports_history_dependence end

supports_history_dependence(::AbstractTimeSolver) = true
supports_history_dependence(::AbstractMaterial) = true

"""
    check_constitutive_model(spatial_setup, solver)

$(internal_api_warning())

Check that a history-dependent constitutive model is combined with a time solver that
evaluates it once per time step, with a material whose evaluation index is well defined,
and with a storage that carries the state of the model. Throws a
[`HistoryDependenceError`](@ref) that names the reason and how to fix it.

This check is done once when a [`Job`](@ref) is created, next to
[`check_storage_contract`](@ref).
"""
function check_constitutive_model(mat::AbstractMaterial, solver::AbstractTimeSolver)
    model = get_constitutive_model(mat)
    is_history_dependent(model) || return nothing
    M, S = typeof(mat), typeof(solver)
    if !supports_history_dependence(solver)
        reason = "the time solver `$(nameof(S))` evaluates the force density more than "
        reason *= "once per time step, so the model would integrate its history several "
        reason *= "times per step"
        fix = "Use a time solver that evaluates the force density once per step, e.g. "
        fix *= "`VelocityVerlet`, or a constitutive model without history."
        throw(HistoryDependenceError(model, M, S, reason, fix))
    end
    if !supports_history_dependence(mat)
        reason = "the material `$(nameof(M))` does not always evaluate the constitutive "
        reason *= "model for the same kind of quantity, so the state of the model has no "
        reason *= "well defined index"
        fix = "Use a material family whose constitutive model is evaluated per bond or per "
        fix *= "point throughout, e.g. `RKCMaterial` or `CMaterial`."
        throw(HistoryDependenceError(model, M, S, reason, fix))
    end
    if constitutive_storage_type(model) !== Nothing &&
       !has_constitutive_state(storage_type(mat))
        reason = "the storage `$(nameof(storage_type(mat)))` of `$(nameof(M))` does not "
        reason *= "carry the state of the constitutive model"
        fix = "Add the state to the `Peridynamics.@storage` definition of "
        fix *= "`$(nameof(storage_type(mat)))`:\n        cm_state::ConstitutiveState"
        throw(HistoryDependenceError(model, M, S, reason, fix))
    end
    return nothing
end

function check_constitutive_model(body::AbstractBody, solver::AbstractTimeSolver)
    return check_constitutive_model(body.mat, solver)
end

function check_constitutive_model(ms::AbstractMultibodySetup, solver::AbstractTimeSolver)
    for body in each_body(ms)
        check_constitutive_model(body, solver)
    end
    return nothing
end
