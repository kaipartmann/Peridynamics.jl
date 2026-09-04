# Extension API

The names on this page are what you need to add a material, a constitutive model, a damage
model or a set of point parameters of your own. They are declared `public` in
`src/public_api.jl`, are not exported and are stable within a minor release.

[Extending Peridynamics.jl](@ref) is the index of the developer documentation, with the rule
the names are built by and the page for each extension point. See [API stability](@ref) for
what the tier promises.

```@meta
CollapsedDocStrings = true
```

```@contents
Pages = ["extension_api_reference.md"]
Depth = 2:2
```

## Declaration macros

The parts of a material are declared, not written by hand.
See [Point parameters](@ref), [Storages](@ref) and [Systems](@ref) for the declaration
language.

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
its element type is and how it is allocated. The table of shapes is in [Storages](@ref).

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

## System types

The discretization of a body chunk, which a material dispatches on. Its fields are internal,
a material reads it through the accessors under
[Accessing a system and its parameters](@ref). See [Systems](@ref).

```@docs
Peridynamics.BondSystem
Peridynamics.InteractionSystem
```

The rest of the system interface, which is internal beyond the two types above:

```@docs
Peridynamics.system_type
Peridynamics.check_system_compat
Peridynamics.SystemSizes
Peridynamics.host_system_type
Peridynamics.max_n_chunks
Peridynamics.first_chunk
```

## The material interface

What a material defines, see [Materials](@ref).

```@docs
Peridynamics.force_density_point!
Peridynamics.storage_type
Peridynamics.init_field
```

## The constitutive model interface

A constitutive model does not depend on the material family that evaluates it, so the same
model runs on `CMaterial`, `RKCMaterial` and `BACMaterial`. See
[Constitutive models](@ref).

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

A damage model decides which bonds fail, and it may carry per-bond state of its own. It is a
plug-in box, and these names are its walls. See [Damage models](@ref).

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

What a force density or a failure criterion reads, one accessor per quantity of a bond. The
lookup table is in [Materials](@ref).

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

Every storage field is one array with the quantity of a point or a bond in its columns, and
these functions read and write a column as a static vector or tensor without allocating. See
[Storages](@ref).

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
Peridynamics.SofteningSupportError
```
