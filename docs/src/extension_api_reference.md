# Extension API

The names on this page are what you need to add a material, a constitutive model, a damage
model or a set of point parameters of your own. They are not exported, so write them as
`Peridynamics.<name>` or import them explicitly.

They are declared `public` in `src/public_api.jl` and are stable within a minor release. See
[API stability](@ref) for what exactly that promises and which parts of the package are
deliberately left out of it. The tutorials [Writing your own material](@ref
tutorial_custom_material), [Writing your own damage model](@ref
tutorial_custom_damage_model) and [Writing your own constitutive model](@ref
tutorial_custom_constitutive_model) show them in use, and [Materials](@ref) is the manual
of the declaration language.

The names follow four rules, so that a name you have not seen yet still tells you what it
does. `get_*` pulls something that is already stored out of a container, as
`get_params` and `get_vector_diff` do. `each_*_idx` returns what a loop iterates over.
`update_*!` writes in place. Everything else is a bare noun, either a physical quantity that
is computed on the spot, as `kernel`, `surface_correction_factor`, `current_bond_length` and
`bond_stretch` are, or a question asked of a material or a model, as `damage_state` and
`storage_type` are.

```@meta
CollapsedDocStrings = true
```

```@contents
Pages = ["extension_api_reference.md"]
Depth = 2:2
```

## Declaration macros

The parts of a material are declared, not written by hand: the point parameters with
[`@params`](@ref Peridynamics.@params), the storage with
[`@storage`](@ref Peridynamics.@storage), the parameters a constitutive or damage model owns
with [`@cm_params`](@ref Peridynamics.@cm_params) and
[`@dmg_params`](@ref Peridynamics.@dmg_params), the state of a model with
[`@cm_storage`](@ref Peridynamics.@cm_storage) and
[`@dmg_storage`](@ref Peridynamics.@dmg_storage), and reusable blocks of declarations with
[`@params_fields`](@ref Peridynamics.@params_fields) and
[`@storage_fields`](@ref Peridynamics.@storage_fields). A system is declared the same way,
with [`@system`](@ref Peridynamics.@system), but has no reusable field blocks of its own: a
system lists every one of its fields directly.

```@docs
Peridynamics.@params
Peridynamics.@params_fields
Peridynamics.@cm_params
Peridynamics.@dmg_params
Peridynamics.@storage
Peridynamics.@storage_fields
Peridynamics.@cm_storage
Peridynamics.@dmg_storage
Peridynamics.@system
```

These are recognized inside the body of the macros above and are never called on their own:

```@docs
Peridynamics.@inherit
Peridynamics.@lth
Peridynamics.@htl
Peridynamics.@kwarg
Peridynamics.@derived
Peridynamics.@log
```

## Field shapes

A storage field is declared with a field shape, which decides how many entries it has, what
its element type is and how it is allocated.

```@docs
Peridynamics.PointScalar
Peridynamics.PointVector
Peridynamics.PointTensor
Peridynamics.PointSymTensor
Peridynamics.PointField
Peridynamics.BondScalar
Peridynamics.BondVector
Peridynamics.BondTensor
Peridynamics.BondSymTensor
Peridynamics.BondField
Peridynamics.DofVector
```

## Field shape modifiers and markers

```@docs
Peridynamics.SimFloat
Peridynamics.LocalPoints
Peridynamics.HaloPoints
Peridynamics.FullField
Peridynamics.EmptyField
Peridynamics.ConstitutiveState
Peridynamics.DamageState
Peridynamics.ConstitutiveParameters
Peridynamics.DamageParameters
```

## Blocks

Everything [`@inherit`](@ref Peridynamics.@inherit) accepts: the parameter and field blocks
of the package, and the point parameters and storages of the shipped materials. What each
of them exposes is generated from its declarations, and the same table is printed by
[`block_table`](@ref Peridynamics.block_table) and by typing the name of a block at the
REPL. [Blocks you can inherit](@ref) is the index of this section.

```@docs
Peridynamics.block_table
```

### Point parameter blocks

The parameter column is what a `@derived` right-hand side can read once the block is
inherited. The keywords are the ones [`material!`](@ref) then accepts.

```@docs
Peridynamics.DiscretizationParameters
Peridynamics.ElasticParameters
Peridynamics.BBElasticParameters
Peridynamics.FractureParameters
Peridynamics.BondHorizonParameters
Peridynamics.InteractionParameters
Peridynamics.StandardParameters
```

### Point parameters of the shipped materials

A point parameter type can be inherited like a block, and it can also be used as it is with
the second form of `@params`, e.g. `Peridynamics.@params MyMaterial BBPointParameters`.

```@docs
Peridynamics.BBPointParameters
Peridynamics.DHBBPointParameters
Peridynamics.OSBPointParameters
Peridynamics.CPointParameters
Peridynamics.RKCPointParameters
Peridynamics.BACPointParameters
Peridynamics.CKIPointParameters
```

### Storage field blocks

A storage needs the block of the time solver it is used with. The fracture bookkeeping
blocks belong to the damage model and are inherited inside a
[`@dmg_storage`](@ref Peridynamics.@dmg_storage) declaration. The others follow from the
system and the material family.

```@docs
Peridynamics.VelocityVerletFields
Peridynamics.DynamicRelaxationFields
Peridynamics.NewtonKrylovFields
Peridynamics.BondLengthCache
Peridynamics.BondFracFields
Peridynamics.InteractionFracFields
Peridynamics.RKCFields
```

### Storages of the shipped materials

A material that builds on a family of this package inherits the storage of that family
instead of listing its fields again, e.g. `@inherit RKCStorage`.

```@docs
Peridynamics.BBStorage
Peridynamics.DHBBStorage
Peridynamics.GBBStorage
Peridynamics.OSBStorage
Peridynamics.CStorage
Peridynamics.CRStorage
Peridynamics.RKCStorage
Peridynamics.RKCRStorage
Peridynamics.BACStorage
Peridynamics.CKIStorage
```

## Abstract types

The supertypes a new material, model or storage is declared with.

```@docs
Peridynamics.AbstractMaterial
Peridynamics.AbstractBondSystemMaterial
Peridynamics.AbstractBondBasedMaterial
Peridynamics.AbstractCorrespondenceMaterial
Peridynamics.AbstractRKCMaterial
Peridynamics.AbstractBondAssociatedSystemMaterial
Peridynamics.AbstractInteractionSystemMaterial
Peridynamics.AbstractConstitutiveModel
Peridynamics.AbstractConstitutiveState
Peridynamics.AbstractDamageModel
Peridynamics.AbstractDamageState
Peridynamics.AbstractStorage
Peridynamics.AbstractPointParameters
Peridynamics.AbstractParameterSetup
Peridynamics.AbstractSystem
Peridynamics.AbstractBondSystem
Peridynamics.AbstractTimeSolver
```

## Systems

The discretization of a body chunk. A material dispatches on it to say which discretization
it is written for. The fields of a system are internal, a material reads it through the
accessors under [Accessing a system and its parameters](@ref), one accessor per quantity of
a bond.

```@docs
Peridynamics.BondSystem
Peridynamics.InteractionSystem
```

Which system a material is discretized with is answered by `system_type`, and
`check_system_compat` is what a system checks its materials with.
`SystemSizes` is what the constructor of a system allocates its fields against, and
`host_system_type` returns the instantiation of a system whose arrays live on the CPU.
A system that cannot be decomposed says so with `max_n_chunks` and reads its one chunk with
`first_chunk`.

```@docs
Peridynamics.system_type
Peridynamics.check_system_compat
Peridynamics.SystemSizes
Peridynamics.host_system_type
Peridynamics.max_n_chunks
Peridynamics.first_chunk
```

## The material interface

```@docs
Peridynamics.force_density_point!
Peridynamics.storage_type
Peridynamics.init_field
```

## The constitutive model interface

A constitutive model does not depend on the material family that evaluates it, so the same
model runs on `CMaterial`, `RKCMaterial` and `BACMaterial`.

```@docs
Peridynamics.get_constitutive_model
Peridynamics.first_piola_kirchhoff
Peridynamics.strain_energy_density
Peridynamics.constitutive_state
Peridynamics.constitutive_storage_type
Peridynamics.is_history_dependent
Peridynamics.supports_history_dependence
```

## The damage model interface

A damage model decides which bonds fail, and it may carry per-bond state of its own. It is
a plug-in box: everything outside of it reads the fracture bookkeeping through
[`bond_is_active`](@ref Peridynamics.bond_is_active) and
[`get_damage`](@ref Peridynamics.get_damage) and never by field name, and a material that
kills bonds writes through [`break_bond!`](@ref Peridynamics.break_bond!) and
[`break_bonds!`](@ref Peridynamics.break_bonds!). The relation between the critical energy
release rate `Gc` and the critical stretch `εc` depends on the micro-modulus, so a material
with a non-constant one defines
[`critical_stretch`](@ref Peridynamics.critical_stretch) and
[`energy_release_rate`](@ref Peridynamics.energy_release_rate). A model that softens a bond
instead of deleting it does so through
[`bond_integrity`](@ref Peridynamics.bond_integrity) for the load the bond still carries and
[`kinematic_weight`](@ref Peridynamics.kinematic_weight) for what it contributes to the
deformation gradient. A material says with
[`supports_bond_integrity`](@ref Peridynamics.supports_bond_integrity) and
[`supports_kinematic_weight`](@ref Peridynamics.supports_kinematic_weight) whether its force
path honors them.

```@docs
Peridynamics.calc_failure!
Peridynamics.calc_damage!
Peridynamics.bond_is_active
Peridynamics.get_damage
Peridynamics.break_bond!
Peridynamics.break_bonds!
Peridynamics.get_dmgmodel
Peridynamics.get_frac_params
Peridynamics.critical_stretch
Peridynamics.energy_release_rate
Peridynamics.has_fracture
Peridynamics.damage_state
Peridynamics.damage_storage_type
Peridynamics.BondFracState
Peridynamics.InteractionFracState
Peridynamics.bond_integrity
Peridynamics.kinematic_weight
Peridynamics.supports_bond_integrity
Peridynamics.supports_kinematic_weight
```

## Accessing a system and its parameters

What a force density or a failure criterion reads: the points and bonds of the chunk
through the iterators, a bond through `get_neighbor`, `reference_bond_length` and
`bond_may_fail`, one accessor per quantity, and the parameters of a point through
`get_params` from the parameter setup the kernel receives.

The current length of a bond and its stretch are read with `current_bond_length` and
`bond_stretch` and never by gathering the two positions and taking the norm. Some materials
keep a cache of the bond lengths and others do not, and these two functions are what makes
the same code as fast as it can be on either.

```@docs
Peridynamics.get_params
Peridynamics.each_point_idx
Peridynamics.each_bond_idx
Peridynamics.get_neighbor
Peridynamics.reference_bond_length
Peridynamics.bond_may_fail
Peridynamics.current_bond_length
Peridynamics.bond_stretch
Peridynamics.update_bond_lengths!
Peridynamics.get_n_points
Peridynamics.get_n_loc_points
Peridynamics.get_n_bonds
Peridynamics.get_n_dim
Peridynamics.kernel
Peridynamics.surface_correction_factor
Peridynamics.float_type
```

## Reading and writing storage fields

Every storage field is one array with the quantity of a point or a bond in its columns.
These functions read and write a column as a static vector or tensor, without allocating.
Each of them takes the number of spatial dimensions as its last argument, a `Val{N}` that
`dims` produces from a system, a storage or a nested model state.

```@docs
Peridynamics.dims
Peridynamics.get_vector
Peridynamics.get_vector_diff
Peridynamics.update_vector!
Peridynamics.update_add_vector!
Peridynamics.get_tensor
Peridynamics.update_tensor!
Peridynamics.get_sym_tensor
Peridynamics.update_sym_tensor!
```

## Strain measures

What a finite strain constitutive model needs in order to work in logarithmic strain space.
Both are evaluated in closed form, without an eigendecomposition and without allocating.

```@docs
Peridynamics.sym_eigvals
Peridynamics.hencky_and_invstretch
```

## Exporting your own fields

```@docs
Peridynamics.export_field
Peridynamics.custom_field
```

## Errors

The errors the interfaces throw when something is missing. They say which method to define.

```@docs
Peridynamics.InterfaceError
Peridynamics.StorageContractError
Peridynamics.HistoryDependenceError
```
