# Constitutive models

The correspondence families do not fix the stress-strain relation. They ask a constitutive
model for the first Piola-Kirchhoff stress that belongs to a deformation gradient, so a new
material behavior usually does not need a new material at all. The tutorial
[Writing your own constitutive model](@ref tutorial_custom_constitutive_model) writes a
hyperelastic model and J2 plasticity in full.

## Contract

| what you define | signature | required | default or fallback |
|:---|:---|:---:|:---|
| the type | `struct MyModel <: Peridynamics.AbstractConstitutiveModel end` | yes | |
| the stress, without a state | [`first_piola_kirchhoff(model, storage, params, F)`](@ref Peridynamics.first_piola_kirchhoff) | one of the two | [`InterfaceError`](@ref Peridynamics.InterfaceError) |
| the stress, with a state | [`first_piola_kirchhoff(model, storage, params, F, idx, Δt)`](@ref Peridynamics.first_piola_kirchhoff) | one of the two | the four-argument form is bridged to this one |
| the parameters of the model | `Peridynamics.@cm_params MyModel struct MyModelParameters ... end` | no | the model has no parameters |
| the state of the model | `Peridynamics.@cm_storage MyModel struct MyModelState ... end` | no | the model has no state and is not history dependent |
| the strain energy density | [`strain_energy_density(model, storage, params, F, idx)`](@ref Peridynamics.strain_energy_density) | only when `:strain_energy_density` is exported | `InterfaceError` |
| that the history is kept elsewhere | [`is_history_dependent(model)`](@ref Peridynamics.is_history_dependent) | no | whether the model declares a state with `@cm_storage` |
| that a solver evaluates the force density once per step | [`supports_history_dependence(solver)`](@ref Peridynamics.supports_history_dependence) | only for a new time solver | `true` |

`strain_energy_density` takes the index but not the time step, and it must not change the
state, because it is also called when a field is exported, that is outside of the time
integration.

## Skeleton

```julia
struct MyModel <: Peridynamics.AbstractConstitutiveModel end

Peridynamics.@cm_params MyModel struct MyModelParameters
    @log "initial yield stress" @kwarg sigma_y σy = Inf
end

Peridynamics.@cm_storage MyModel struct MyModelState
    bond_plastic_strain::BondSymTensor
    bond_eqps::BondScalar
end

function Peridynamics.first_piola_kirchhoff(::MyModel, storage, params, F, idx, Δt)
    state = constitutive_state(storage)
    εp = get_sym_tensor(state.bond_plastic_strain, idx, dims(storage))
    ...
end
```

## What `idx` indexes

| material family | `idx` | state shapes |
|:---|:---|:---|
| [`CMaterial`](@ref) | point index | `Point...` |
| [`RKCMaterial`](@ref), [`BACMaterial`](@ref) | bond index | `Bond...` |

## Which materials carry a state

| material | declares `cm_state::ConstitutiveState` | history-dependent model |
|:---|:---:|:---:|
| [`CMaterial`](@ref) | yes | yes |
| [`RKCMaterial`](@ref) | yes | yes |
| [`BACMaterial`](@ref) | yes | yes |
| [`CRMaterial`](@ref) | no | no |
| [`RKCRMaterial`](@ref) | no | no |

## Rules

- Declaring a state with `@cm_storage` makes the model history dependent, see
  [`is_history_dependent`](@ref Peridynamics.is_history_dependent).
- A history-dependent model may only be run by a time solver that evaluates the force
  density **once** per time step. [`NewtonKrylov`](@ref) evaluates it several times, for the
  Jacobian-vector products and the line search.
- The state is chunk-local and is never exchanged between chunks, which is why the halo
  annotations are not allowed in a `@cm_storage` definition.
- The parameters of a model become keywords of [`material!`](@ref) for every material whose
  point parameters carry the marker `cm_params::ConstitutiveParameters`, and they are read
  flat off the point parameters, e.g. `params.sigma_y`. A declaration may read every
  material parameter declared above the marker, e.g. the shear modulus `μ`, and the model
  instance is available as `model` and the material as `mat`.
- Both the solver and the storage are checked once when a [`Job`](@ref) is created, and a
  mismatch throws a
  [`HistoryDependenceError`](@ref Peridynamics.HistoryDependenceError) that names the
  reason.

See also [Materials](@ref), [Storages](@ref), [Point parameters](@ref),
[Time solvers](@ref), [Blocks you can inherit](@ref), [Extension API](@ref).
