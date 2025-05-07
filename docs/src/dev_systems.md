# Systems

!!! warning "Work in progress"
    This page is currently worked on.

Systems are the backbone of the simulations.
They contain all the required pre-defined data structures (e.g. bonds, point-families, interactions, ...).

## BondSystem

The standard system for bond-based ([`BBMaterial`](@ref)), ordinary state-based ([`OSBMaterial`](@ref)) and non-ordinary state-based correspondence formulation ([`CMaterial`](@ref)) material models. Used material models have to be subtypes of `AbstractBondSystemMaterial` or `AbstractCorrespondenceMaterial`.

## InteractionSystem

The system utilizing one-, two- and three-neighbor-interactions for the continuum-kinematics-inspired material model ([`CKIMaterial`](@ref)).
Used material models have to be subtypes of  `AbstractInteractionSystemMaterial`.

## BondAssociatedSystem

A system for the bond-associated correspondence model by Chen and Spencer ([`BACMaterial`](@ref)).
Used material models have to be subtypes of  `AbstractBondAssociatedSystemMaterial`.

## Custom systems

Custom systems are relatively free in how they are defined, only these things are required:
- define an abstract type `Abstract<SystemName>Material` that all materials using this system have to be a subtype of.
- define the `get_system(body::AbstractBody{Material}, pd::PointDecomposition, chunk_id::Int) where {Material<:Abstract<SystemName>Material}` function.
- define the `system_type(mat::Abstract<SystemName>Material)` function that returns the system type.
- define the `calc_timestep_point(system::<SystemType>, params::AbstractPointParameters, point_id::Int)` function.
- define the `calc_force_density!(chunk::AbstractBodyChunk{<:<SystemType>}, t, Δt)` function.
