# Materials

A material decides which peridynamic formulation a body is simulated with. The tutorial
[Writing your own material](@ref tutorial_custom_material) writes one in full, this page is
the lookup table.

## The materials of this package

| material | formulation | system | constitutive model |
|:---|:---|:---|:---:|
| [`BBMaterial`](@ref) | bond-based | `BondSystem` | no |
| [`DHBBMaterial`](@ref) | dual-horizon bond-based | `BondSystem` | no |
| [`GBBMaterial`](@ref) | generalized bond-based | `BondSystem` | no |
| [`OSBMaterial`](@ref) | ordinary state-based, the linear peridynamic solid | `BondSystem` | no |
| [`CMaterial`](@ref) | correspondence | `BondSystem` | yes |
| [`CRMaterial`](@ref) | correspondence with stress rotation | `BondSystem` | yes |
| [`RKCMaterial`](@ref) | reproducing kernel with bond-associated integration | `BondSystem` | yes |
| [`RKCRMaterial`](@ref) | reproducing kernel with stress rotation | `BondSystem` | yes |
| [`BACMaterial`](@ref) | bond-associated correspondence of Chen and Spencer | `BondAssociatedSystem` | yes |
| [`CKIMaterial`](@ref) | continuum-kinematics-inspired | `InteractionSystem` | no |

## Contract

| what you define | signature | required | default or fallback |
|:---|:---|:---:|:---|
| the supertype, which picks the system | `struct MyMaterial <: Peridynamics.AbstractBondSystemMaterial{Correction}` | yes | none, see the table of supertypes below |
| the point parameters | `Peridynamics.@params MyMaterial struct MyPointParameters ... end` or `Peridynamics.@params MyMaterial BBPointParameters` | yes | [`InterfaceError`](@ref Peridynamics.InterfaceError) |
| the storage | `Peridynamics.@storage MyMaterial struct MyStorage ... end` | yes | [`InterfaceError`](@ref Peridynamics.InterfaceError) |
| the force density | [`force_density_point!(storage, system, mat, paramsetup, t, Δt, i)`](@ref Peridynamics.force_density_point!) | yes | none, a `MethodError` |
| the point parameter `bc` | `@derived bc = ...` inside the `@params` body | under [`VelocityVerlet`](@ref) | `@inherit StandardParameters` derives it, see the rules below |
| a field named `dmgmodel` | `struct MyMaterial{Correction,DM} ... dmgmodel::DM ... end` | only to run a damage model | [`get_dmgmodel`](@ref Peridynamics.get_dmgmodel) returns `nothing` and the material has no fracture |
| an unshaped storage field | [`init_field(mat, solver, system, ::Val{:field})`](@ref Peridynamics.init_field) | only for a field no shape describes | a shaped field is allocated by its shape |
| a field of your own in the VTK output | [`custom_field(::Type{<:MyStorage}, ::Val{:name})`](@ref Peridynamics.custom_field) and [`export_field(::Val{:name}, mat, system, storage, paramsetup, t)`](@ref Peridynamics.export_field) | only for a derived quantity | every point field of the storage is exported as it is |
| another micro-modulus | [`critical_stretch(dmgmodel, mat, δ, K, Gc)`](@ref Peridynamics.critical_stretch) and [`energy_release_rate(dmgmodel, mat, δ, K, εc)`](@ref Peridynamics.energy_release_rate), always both | only with a non-constant micro-modulus | the relation of the constant micro-modulus |
| that the force path softens bonds | [`supports_bond_integrity(mat)`](@ref Peridynamics.supports_bond_integrity), [`supports_kinematic_weight(mat)`](@ref Peridynamics.supports_kinematic_weight) | no | `false` |
| that the system cannot be decomposed | [`max_n_chunks(mat)`](@ref Peridynamics.max_n_chunks) | no | `typemax(Int)` |

The supertype is what says which system a material is discretized on:

| supertype | system |
|:---|:---|
| [`AbstractBondSystemMaterial{Correction}`](@ref Peridynamics.AbstractBondSystemMaterial) | `BondSystem` |
| [`AbstractCorrespondenceMaterial`](@ref Peridynamics.AbstractCorrespondenceMaterial) | `BondSystem`, with a constitutive model |
| [`AbstractRKCMaterial`](@ref Peridynamics.AbstractRKCMaterial) | `BondSystem`, reproducing kernel |
| [`AbstractBondAssociatedSystemMaterial`](@ref Peridynamics.AbstractBondAssociatedSystemMaterial) | `BondAssociatedSystem` |
| [`AbstractInteractionSystemMaterial`](@ref Peridynamics.AbstractInteractionSystemMaterial) | `InteractionSystem` |

## Skeleton

```julia
struct MyMaterial{Correction,DM} <: Peridynamics.AbstractBondSystemMaterial{Correction}
    dmgmodel::DM
end

function MyMaterial{C}(; dmgmodel=CriticalStretch()) where {C}
    return MyMaterial{C,typeof(dmgmodel)}(dmgmodel)
end
MyMaterial(; kwargs...) = MyMaterial{NoCorrection}(; kwargs...)

Peridynamics.@params MyMaterial struct MyPointParameters
    @inherit StandardParameters
    ...
end

Peridynamics.@storage MyMaterial struct MyStorage
    @inherit VelocityVerletFields DynamicRelaxationFields NewtonKrylovFields
    @inherit BondLengthCache
    dmg_state::DamageState
    ...
end

function Peridynamics.force_density_point!(storage::MyStorage, system::BondSystem,
                                           mat::MyMaterial, paramsetup, t, Δt, i)
    ...
end
```

## What you may read

| function | returns | notes |
|:---|:---|:---|
| [`each_bond_idx(system, i)`](@ref Peridynamics.each_bond_idx) | the bond indices of point `i` | addresses every bond field of the system and the storage |
| [`get_neighbor(system, bond_id)`](@ref Peridynamics.get_neighbor) | the point index `j` of the bond | |
| [`reference_bond_length(system, bond_id)`](@ref Peridynamics.reference_bond_length) | the initial length `L` | |
| [`bond_may_fail(system, bond_id)`](@ref Peridynamics.bond_may_fail) | whether the bond is allowed to break | how [`no_failure!`](@ref) is honored |
| [`kernel(system, bond_id)`](@ref Peridynamics.kernel) | the influence function of the bond | |
| [`surface_correction_factor(system, bond_id)`](@ref Peridynamics.surface_correction_factor) | the correction factor | `1` with `NoCorrection` |
| `system.volume[j]` | the volume of the neighbor | the one system field read by name |
| [`get_vector_diff(storage.position, i, j, dims(system))`](@ref Peridynamics.get_vector_diff) | the bond vector as an `SVector{N}` | works on `system.position` for the reference vector |
| [`current_bond_length(storage, system, i, bond_id)`](@ref Peridynamics.current_bond_length) | the deformed length `l` | reads the cache where there is one |
| [`bond_stretch(storage, system, i, bond_id)`](@ref Peridynamics.bond_stretch) | the stretch `ε` | one read, not length and reference length |
| [`bond_is_active(storage, system, bond_id)`](@ref Peridynamics.bond_is_active) | whether the bond is intact | the damage model decided this right before |
| [`get_params(paramsetup, i)`](@ref Peridynamics.get_params) | the point parameters of point `i` | free on a body with one parameter set |
| [`update_add_vector!(storage.b_int, i, b, dims(system))`](@ref Peridynamics.update_add_vector!) | writes, adds `b` to column `i` | how a force density accumulates |

The whole force density of a bond-based material reads like this:

```julia
function Peridynamics.force_density_point!(storage::MyStorage, system::BondSystem,
                                           mat::MyMaterial, paramsetup, t, Δt, i)
    (; volume) = system
    params = get_params(paramsetup, i)
    for bond_id in each_bond_idx(system, i)
        j = get_neighbor(system, bond_id)
        L = reference_bond_length(system, bond_id)
        Δxij = get_vector_diff(storage.position, i, j, dims(system))
        l = current_bond_length(storage, system, i, bond_id)
        ε = (l - L) / L
        ω = bond_is_active(storage, system, bond_id) *
            surface_correction_factor(system, bond_id)
        b = ω * params.bc * ε * volume[j] / l .* Δxij
        update_add_vector!(storage.b_int, i, b, dims(system))
    end
    return nothing
end
```

## Rules

- The kernel writes **only columns of point `i`**. Everything else is read.
- Read the deformed length with `current_bond_length` and the stretch with `bond_stretch`,
  never by gathering the two positions and taking the norm. Some materials cache the bond
  lengths and others do not, and these two functions are what makes the same line as fast
  as it can be either way. A kernel that needs the length forms the stretch from it, so
  that the bond is read once.
- Ask `bond_is_active` instead of reading a storage field. That is what makes the same
  kernel run with every damage model, including one that carries no bookkeeping.
- A body may have several parameter sets, one per point set. `get_params(paramsetup, i)`
  resolves the set of point `i` either way. A material that averages a parameter over the
  two points of a bond reads `get_params(paramsetup, j)` inside the loop. On a body with a
  single set that read does not depend on the loop, so the averaging moves out of it and
  costs nothing. `BBMaterial`, `DHBBMaterial`, `GBBMaterial`, `OSBMaterial` and
  `CKIMaterial` are all written with this one kernel.
- **The name `bc` is not free.** The stable time step of an explicit solver is estimated
  from the bond constant, so a material that derives its own has to keep the name. If the
  bond stiffness is not constant over the family, declare `bc` as its largest value, so
  that the estimate stays on the safe side. Such a material also defines `critical_stretch`
  and `energy_release_rate`, because the relation between `Gc` and `εc` follows from the
  micro-modulus.
- **The name `dmgmodel` is not free** either. Every bond system material is asked for its
  damage model before the force density is evaluated, and `get_dmgmodel` looks for that
  field by name.
- The complete storage contract is checked once when a [`Job`](@ref) is created, and a
  missing field throws a
  [`StorageContractError`](@ref Peridynamics.StorageContractError) that names the field and
  the reason why it is required.

## Exporting a field of your own

Every point field of a storage can be named in the `fields` keyword of a [`Job`](@ref) and
is written to the VTK files as it is. A quantity that is not a storage field, or a bond
field that has to be reduced to one value per point, is announced with `custom_field` and
computed by `export_field`:

```julia
Peridynamics.custom_field(::Type{<:MyStorage}, ::Val{:bond_damage_avg}) = true

function Peridynamics.export_field(::Val{:bond_damage_avg}, mat, system, storage::MyStorage,
                                   paramsetup, t)
    ...
end
```

The returned array has one entry per **local** point, i.e.
`Peridynamics.get_n_loc_points(system)` of them. Halo entries belong to another chunk and
must not be exported twice.

See also [Point parameters](@ref), [Storages](@ref), [Damage models](@ref),
[Constitutive models](@ref), [Systems](@ref), [Blocks you can inherit](@ref),
[Extension API](@ref).
