# Materials

!!! warning "Draft & Work in Progress"
    This documentation is a draft and work in progress. It will be extended and improved in the future.

The selected material specifies, which peridynamic formulation is employed in calculations.

## Implemented material models

- [`BBMaterial`](@ref): Bond-based peridynamics.
- [`DHBBMaterial`](@ref): Dual-horizon bond-based peridynamics.
- [`OSBMaterial`](@ref): Ordinary state-based peridynamics, also called linear peridynamic solid (LPS).
- [`CMaterial`](@ref): Correspondence formulation.
- [`CRMaterial`](@ref): Correspondence formulation with stress rotation for objectivity enforcement.
- [`RKCMaterial`](@ref): Reproducing kernel peridynamics with bond-associated higher-order integration.
- [`RKCRMaterial`](@ref): Reproducing kernel peridynamics with bond-associated higher-order integration with stress rotation for objectivity enforcement.
- [`BACMaterial`](@ref): Bond-associated correspondence formulation of Chen and Spencer.
- [`CKIMaterial`](@ref): Continuum-kinematics-inspired peridynamics.

## Custom materials

Custom materials can be defined for existing systems.
Therefore the `<XYZ>Material` type has to be a subtype of an `Abstract<SystemName>Material`, which automatically sets methods for this type.
Further it is necessary to: 
- link a point parameter type with `@params`.
- link and create a storage with `@storage`.
- define the `force_density_point!` function.
