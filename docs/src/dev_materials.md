# Materials

!!! warning "Work in progress"
    This page is currently worked on.

The selected material specifies, which peridynamic formulation is employed in calculations.

## Implemented material models

- [`BBMaterial`](@ref): Bond-based material model.
- [`OSBMaterial`](@ref): Ordinary state-based material model.
- [`CMaterial`](@ref): Non-ordinary state-based corespondence material model.
- [`CKIMaterial`](@ref): Continuum-kinematics-inspired material model.
- [`BACMaterial`](@ref): Bond-associated correspondence material model by Chen and Spencer.

## Custom materials

Custom materials can be defined for existing systems.
Therefore the `<XYZ>Material` type has to be a subtype of an `Abstract<SystemName>Material`, which automatically sets methods for this type.
Further it is necessary to: 
- link a point parameter type with `@params`.
- link and create a storage with `@storage`.
- define the `force_density_point!` function.
