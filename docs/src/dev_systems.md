# Systems

!!! warning "Work in progress"
    This page is currently worked on.

## AbstractBondSystem

## AbstractInteractionSystem

## Custom systems

Systems are the backbone of the simulations
They contain all the required pre-defined data structures (e.g. bonds, point-families, interactions, ...)
BondSystem: The standard system for BB-, OSB-, C-Material
InteractionSystem: The system utilizing one-, two-, three-neighbor-interactions for the CKI-Material
They are relatively free in how they are defined, only these things are required:
define a abstract type Abstract<SystemName>Material that all material using this system have to be a subtype of
define the get_system(body::AbstractBody{Material}, pd::PointDecomposition, chunk_id::Int) where {Material<:Abstract<SystemName>Material} function
define the system_type(mat::Abstract<SystemName>Material) function that returns the system type
define the calc_timestep_point(system::<SystemType>, params::AbstractPointParameters, point_id::Int) function
define the calc_force_density!(chunk::AbstractBodyChunk{<:<SystemType>}, t, Δt) function
... Work in progress ...