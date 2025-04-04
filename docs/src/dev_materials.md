# Materials

!!! warning "Work in progress"
    This page is currently worked on.

## Implemented material models

- [`BBMaterial`](@ref)
- [`OSBMaterial`](@ref)
- [`CMaterial`](@ref)
- [`CKIMaterial`](@ref)
- [`BACMaterial`](@ref)

## Custom materials

Custom materials can be defined for existing systems.
Therefore the `<XYZ>Material` type has to be subtype of an `Abstract<Systemname>Material`, which automatically sets methods for this type.
Further it is necessary to: 
- link a point parameter type with `@params`.
- link and create a storage with `@storage`.
- define the `force_density_point!` function.

