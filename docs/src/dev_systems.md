# Systems

A system is the discretization of a body chunk: its points, their volumes and their
neighborhood relations. It is built once during setup and never changes during a
simulation, which is what separates it from the storage.

!!! warning "Writing a new system is an internal interface"
    [`BondSystem`](@ref Peridynamics.BondSystem) and
    [`InteractionSystem`](@ref Peridynamics.InteractionSystem) are part of the
    [Extension API](@ref), because a material has to dispatch on them. The interface for
    adding a new system is not. The names in the last section can change in any release.
    It will be stabilized in a later release, so if you are writing one, please open an
    issue.

## The systems of this package

### BondSystem

The standard system, in which the neighborhood relation is a bond between two points. It is
used by the bond-based ([`BBMaterial`](@ref), [`DHBBMaterial`](@ref), [`GBBMaterial`](@ref)),
the ordinary state-based ([`OSBMaterial`](@ref)), the correspondence ([`CMaterial`](@ref),
[`CRMaterial`](@ref)) and the reproducing kernel ([`RKCMaterial`](@ref),
[`RKCRMaterial`](@ref)) materials. Materials are subtypes of
[`AbstractBondSystemMaterial`](@ref Peridynamics.AbstractBondSystemMaterial).

Inside a force density calculation the points and bonds of a chunk are walked with the
iterators, and a bond is read off `system.bonds`:

```julia
for i in each_point_idx(system)                        # the local points of this chunk
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length              # neighbor index, initial length
        ωij = kernel(system, bond_id)                  # the influence function
        β = surface_correction_factor(system.correction, bond_id)  # 1 with `NoCorrection`
        Vj = system.volume[j]                          # the volume of the neighbor
        ΔXij = get_vector_diff(system.position, i, j)  # the initial bond vector
    end
end
```

A [`Bond`](@ref Peridynamics.Bond) also carries `fail_permit`, which a damage model reads.
[`get_n_loc_points`](@ref Peridynamics.get_n_loc_points) is the number of points this chunk
integrates, [`get_n_points`](@ref Peridynamics.get_n_points) additionally counts the halo
points that are read from neighboring chunks, and
[`get_n_bonds`](@ref Peridynamics.get_n_bonds) is the number of bonds.

### BondAssociatedSystem

A bond system that additionally carries the bond-associated neighborhoods used by the
correspondence model of Chen and Spencer ([`BACMaterial`](@ref)). Materials are subtypes of
[`AbstractBondAssociatedSystemMaterial`](@ref Peridynamics.AbstractBondAssociatedSystemMaterial).

### InteractionSystem

A system built on one-, two- and three-neighbor interactions instead of bonds, used by the
continuum-kinematics-inspired material ([`CKIMaterial`](@ref)). Materials are subtypes of
[`AbstractInteractionSystemMaterial`](@ref Peridynamics.AbstractInteractionSystemMaterial).

## Local points, halo points and the field extents

A body is decomposed into chunks, one per thread or per MPI rank. Each chunk owns its
**local points** and integrates their equation of motion. To do that it needs the current
state of some points of other chunks, its **halo points**, which are exchanged before and
after the force density calculation.

This is the whole reason a storage field declares an extent:

- a field with the default extent [`LocalPoints`](@ref Peridynamics.LocalPoints) has
  `get_n_loc_points(system)` entries,
- a field annotated with [`@lth`](@ref Peridynamics.@lth) or
  [`@htl`](@ref Peridynamics.@htl) has `get_n_points(system)` entries, so it has halo
  entries and is exchanged.

`@lth` copies the local entries of the owning chunk into the halo entries of its neighbors,
which is what `position` needs. `@htl` adds the halo entries back into the local entries of
the owner, which is what a material that accumulates a force density into its neighbors
needs (`@htl b_int::PointVector`).

A material never has to do anything for this beyond the annotation. The exchange is the
same code under multithreading and under MPI.

## Custom systems

A custom system is relatively free in how it is defined. What is required:

- an abstract type `AbstractMySystemMaterial` that every material using this system is a
  subtype of,
- `Peridynamics.get_system(body::AbstractBody{Material}, pd::PointDecomposition, chunk_id::Int)`
  for `Material <: AbstractMySystemMaterial`,
- `Peridynamics.system_type(mat::AbstractMySystemMaterial)`, returning the system type,
- `Peridynamics.calc_timestep_point(system::MySystem, params, point_id::Int)`,
- `Peridynamics.calc_force_density!(chunk::AbstractBodyChunk{<:MySystem}, t, Δt)`,
- `Peridynamics.float_type(::MySystem)` and the chunk handler accessors, so that the
  extents above resolve.
