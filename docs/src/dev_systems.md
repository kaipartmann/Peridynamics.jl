# Systems

A system is the discretization of a body chunk: its points, their volumes and their
neighborhood relations. It is built once during setup and never changes during a
simulation, which is what separates it from the storage.

!!! warning "Writing a new system is an internal interface"
    [`BondSystem`](@ref Peridynamics.BondSystem) and
    [`InteractionSystem`](@ref Peridynamics.InteractionSystem) are part of the
    [Extension API](@ref), because a material has to dispatch on them. The interface for
    adding a new system is not. The names in the last sections can change in any release.
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

### BondAssociatedSystem

A bond system that additionally carries the bond-associated neighborhoods used by the
correspondence model of Chen and Spencer ([`BACMaterial`](@ref)). Materials are subtypes of
[`AbstractBondAssociatedSystemMaterial`](@ref Peridynamics.AbstractBondAssociatedSystemMaterial).
The family of a bond is walked with `each_intersecting_bond_idx(system, i, bond_id)`, which
gives the bond indices of the chunk, so the same bond fields are addressed inside and
outside the family.

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

* a field with the default extent [`LocalPoints`](@ref Peridynamics.LocalPoints) has
  `get_n_loc_points(system)` entries,
* a field annotated with [`@lth`](@ref Peridynamics.@lth) or
  [`@htl`](@ref Peridynamics.@htl) has `get_n_points(system)` entries, so it has halo
  entries and is exchanged.

`@lth` copies the local entries of the owning chunk into the halo entries of its neighbors,
which is what `position` needs. `@htl` adds the halo entries back into the local entries of
the owner, which is what a material that accumulates a force density into its neighbors
needs (`@htl b_int::PointVector`).

A material never has to do anything for this beyond the annotation. The exchange is the
same code under multithreading and under MPI.

## Reading a system in a kernel

Whatever a system is built from, a material reads it through functions and never by field
name. The points of a chunk are walked with
[`each_point_idx`](@ref Peridynamics.each_point_idx), the bonds of a point with
[`each_bond_idx`](@ref Peridynamics.each_bond_idx), and a bond is read one quantity at a
time, each of them a single array load:

```julia
for i in each_point_idx(system)                         # the local points of this chunk
    for bond_id in each_bond_idx(system, i)
        j = get_neighbor(system, bond_id)               # the neighbor index
        L = reference_bond_length(system, bond_id)      # the initial length
        ωij = kernel(system, bond_id)                   # the influence function
        β = surface_correction_factor(system, bond_id)  # 1 with `NoCorrection`
        Vj = system.volume[j]                           # the volume of the neighbor
        ΔXij = get_vector_diff(system.position, i, j, dims(system))  # the bond vector
    end
end
```

[`bond_may_fail`](@ref Peridynamics.bond_may_fail) is what a damage model asks instead of
reading a `fail_permit` field, and
[`current_bond_length`](@ref Peridynamics.current_bond_length) and
[`bond_stretch`](@ref Peridynamics.bond_stretch) are how the deformed length is read, so
that a system with a length cache and one without behave the same.

The last argument of every accessor of a storage or system field is the number of spatial
dimensions as a `Val`, which [`dims`](@ref Peridynamics.dims) produces from whatever is in
scope: `dims(system)`, `dims(storage)` or `dims(state)`. Both
[`get_n_dim`](@ref Peridynamics.get_n_dim) and [`float_type`](@ref Peridynamics.float_type)
read a type parameter of the system, so the `Val` is a compile time constant and the
accessor folds into plain indexing. `dims(body)` is not a constant and must never appear in
a kernel.

## The system contract

This is everything a system is asked for. A system generated by
[`@system`](@ref Peridynamics.@system) gets the rows marked "generated" for free, the rest
is written by hand and is one line each in most cases.

| function | who calls it | required | generated |
|:---|:---|:---:|:---:|
| [`system_type(mat, FT, Val(N))`](@ref Peridynamics.system_type) | `body_chunk_type`, `damage_storage_type`, `log_system` | yes | no |
| `get_system(body, pd, chunk_id)` and the constructor | `BodyChunk` | yes | no |
| [`get_n_dim`](@ref Peridynamics.get_n_dim), [`float_type`](@ref Peridynamics.float_type) | the shape layer, the dof helpers | yes | yes |
| [`get_n_loc_points`](@ref Peridynamics.get_n_loc_points), [`get_n_points`](@ref Peridynamics.get_n_points), [`each_point_idx`](@ref Peridynamics.each_point_idx), `get_point_ids`, `get_localizer`, `get_hidxs_by_src`, `get_loc_view` | the shape layer, the halo exchange, pre-cracks, conditions, export | yes | through `chunk_handler` |
| [`get_n_bonds`](@ref Peridynamics.get_n_bonds) | bond shaped storage fields | if it has bonds | yes |
| `calc_force_density!(chunk, t, Δt)` | the data handlers | yes | the family level exists |
| `calc_timestep_point(system, params, i)` | [`VelocityVerlet`](@ref) | yes | no |
| [`calc_damage!(chunk)`](@ref Peridynamics.calc_damage!) | `apply_precracks!` | yes | the family level exists |
| `log_system(::Type{<:MySystem}, options, dh)` | the data handlers | yes | no |
| `Adapt.adapt_structure`, [`host_system_type`](@ref Peridynamics.host_system_type) | [`system_type`](@ref Peridynamics.system_type), the device | yes | yes |
| `init_field_system(system, Val(f))` | storage construction | optional, `position` is generic | no |
| `initialize!(chunk)`, `initialize!(dh, solver)` | the data handlers | optional | no |
| `apply_precrack!(chunk, body, crack)` | `apply_precracks!` | optional | no |
| [`max_n_chunks(mat)`](@ref Peridynamics.max_n_chunks) | `threads_data_handler`, the MPI check | optional, default `typemax(Int)` | no |
| [`first_chunk(dh)`](@ref Peridynamics.first_chunk) | `log_system` of a system that exists once per body | the data handlers provide it | not applicable |
| [`each_bond_idx`](@ref Peridynamics.each_bond_idx), [`get_neighbor`](@ref Peridynamics.get_neighbor), [`reference_bond_length`](@ref Peridynamics.reference_bond_length), [`bond_may_fail`](@ref Peridynamics.bond_may_fail), [`kernel`](@ref Peridynamics.kernel), [`surface_correction_factor`](@ref Peridynamics.surface_correction_factor), [`current_bond_length`](@ref Peridynamics.current_bond_length), [`bond_stretch`](@ref Peridynamics.bond_stretch), [`update_bond_lengths!`](@ref Peridynamics.update_bond_lengths!) | kernels, damage models | if it has bonds | shared through the contract of [`AbstractBondSystem`](@ref Peridynamics.AbstractBondSystem) |
| [`check_system_compat(::Type{MySystem}, mat)`](@ref Peridynamics.check_system_compat) | the constructors | optional | no |

## Declaring a system

A system is declared with [`@system`](@ref Peridynamics.@system), which is
[`@storage`](@ref Peridynamics.@storage) for a system. The body accepts the same field
shapes, so a field of a system is sized, allocated and moved to another array backend
exactly like a storage field:

```julia
Peridynamics.@system struct MySystem <: Peridynamics.AbstractSystem
    position::PointVector{Float64}
    volume::PointScalar
    neighbor::BondScalar{Int}
    bond_length::BondScalar
    bond_ids::PointScalar{UnitRange{Int}}
end
```

A system has no reusable field blocks: `@storage_fields` and `@inherit` are for a storage,
and writing `@inherit` inside `@system` is an error. A system lists every one of its fields
directly, so the three systems of this package (`BondSystem`, `BondAssociatedSystem`,
`InteractionSystem`) each repeat the bond fields `position`, `volume`, `neighbor`,
`bond_length`, `fail_permit`, `n_neighbors` and `bond_ids` at the top of their declaration.

Two rules decide what a declaration becomes:

* a field declared with a **field shape** gets the container, the element type and the
  number of entries of that shape, so `volume::PointScalar` is a vector of the float type
  with one entry per point and `neighbor::BondScalar{Int}` is a vector of `Int` with one
  entry per bond,
* every other declared type is kept as it is. An isbits scalar or a struct of the system,
  e.g. `lattice::FastLattice`, is left alone by `Adapt` and is not a parameter unless it is
  a concrete `Array` type, which is treated like a shaped field. A field declared with an
  **abstract type** is an error: give it a concrete type, or declare a type parameter of
  the system bounded by that abstract type and use that for the field.

Unlike a storage field, a system field never declares an initial value with `= value`. The
constructor of the system fills every field, so an initial value would be written and never
read.

A system never declares a `chunk_handler` field itself, doing so is an error. The macro
injects it as the last field of the struct, which is what answers the point counts and the
halo bookkeeping of the section above.

### What the macro generates

The type parameters of the generated struct are the ones the system declares itself, then
`N` for the number of spatial dimensions, then `FT` for the float type of the simulation,
then `CH` for the type of the injected chunk handler, then one parameter per distinct array
type of its fields. The declared parameters come first so that a pattern like
`BondSystem{<:EnergySurfaceCorrection}` keeps selecting a method after `N`, `FT`, `CH` and
the array parameters were appended behind it. A user type parameter named `N`, `FT` or `CH`
is an error, because the macro fills those in itself.

The macro generates

* the positional constructor `MySystem{N,FT}(fields...,chunk_handler)`, which takes the
  declared fields in the order they were declared and the chunk handler last, and infers
  `CH` and the array parameters from the values given, and checks that `position` really
  has `N` rows,
* `Adapt.adapt_structure`,
* [`host_system_type`](@ref Peridynamics.host_system_type), the instantiation whose arrays
  live on the CPU,
* [`get_n_dim`](@ref Peridynamics.get_n_dim) and
  [`float_type`](@ref Peridynamics.float_type), which read `N` and `FT` off the type,
* [`get_n_bonds`](@ref Peridynamics.get_n_bonds), but only for a system that declares a
  bond shaped field, so a system without bonds never gets one.

The system is the authority on `N` and `FT`. A storage is built with the `N` of its system,
so the two can never disagree.

### Allocating the fields

A constructor allocates against [`SystemSizes`](@ref Peridynamics.SystemSizes), which is
what a system looks like before it exists. It answers the questions the shape layer asks,
the dimension, the float type, the point counts and the number of bonds, and nothing else,
so `alloc_field` sizes a system field exactly as it sizes a storage field:

```julia
sizes = Peridynamics.SystemSizes{N,FT}(chunk_handler, length(neighbor))
kernels = Peridynamics.alloc_field(BondScalar(), sizes, Peridynamics.LocalPoints())
```

This is what replaces a hand written `zeros(3, n_points)`, and it is why a shaped field of a
system follows the float type and the dimension without the constructor naming either.

## Dispatching on a system

Every dispatch on a system has to be written with `<:`, e.g.
`AbstractBodyChunk{<:MySystem}` and `{<:BondSystem{<:EnergySurfaceCorrection}}`. A pattern
that names the parameters exactly, e.g. `AbstractBodyChunk{MySystem}`, stops matching as
soon as a parameter is added, and it does so silently, because the fallback method it then
selects is type stable too. Trailing parameters may be left out of a `<:` pattern, which is
what makes `BondSystem{<:EnergySurfaceCorrection}` readable.

## Moving a chunk to another array backend

Every array of a system is behind a type parameter, so a whole body chunk moves with
`Adapt.adapt`:

```julia
device_chunk = Adapt.adapt(CuArray, chunk)
```

The system, the parameter setup and the storage move. The material, the conditions and the
export cells stay as they are, because they are host data that no kernel reads. The
`ChunkHandler` becomes a `DeviceChunkHandler`, which keeps the two point counts and leaves
the point ids, the halo bookkeeping and the localizer on the host where the halo exchange
and the export read them. `get_n_dim` and `float_type` of the moved system answer exactly
what they answered before, and a target that moves no array at all, e.g.
`Adapt.adapt(Array, chunk)` on the host, gives the chunk back unchanged.

## A system that cannot be decomposed

Some systems exist once per body, for example one that transforms the whole body at once
and therefore needs every point of it in one place. Such a system says so with
[`max_n_chunks`](@ref Peridynamics.max_n_chunks) on its material:

```julia
Peridynamics.max_n_chunks(::MyMaterial) = 1
```

A threaded run clamps the number of chunks to this value, so nothing has to be configured
by the user. An MPI run throws instead, because the number of ranks comes from the outside
and cannot be clamped away. Where such a system needs its one chunk, for example in
`log_system`, it reads it with [`first_chunk`](@ref Peridynamics.first_chunk), which both
data handlers provide.

## What is left to write

Beyond the declaration, a new system needs

* an abstract type `AbstractMySystemMaterial` that every material using this system is a
  subtype of,
* the constructor `MySystem(body, pd, chunk_id)`, which allocates its fields against a
  [`SystemSizes`](@ref Peridynamics.SystemSizes) as shown above,
* `Peridynamics.get_system(body::AbstractBody{Material}, pd::PointDecomposition, chunk_id::Int)`
  for `Material <: AbstractMySystemMaterial`,
* [`Peridynamics.system_type`](@ref Peridynamics.system_type), one line over
  [`host_system_type`](@ref Peridynamics.host_system_type):

  ```julia
  function Peridynamics.system_type(mat::AbstractMySystemMaterial,
                                    ::Type{FT}=Peridynamics.default_float_type(),
                                    ::Val{N}=Val(3)) where {FT,N}
      return Peridynamics.host_system_type(MySystem, Val(N), FT)
  end
  ```

* [`Peridynamics.check_system_compat`](@ref Peridynamics.check_system_compat), if the system
  accepts only one family of materials,
* `Peridynamics.calc_timestep_point(system::MySystem, params, point_id::Int)`,
* `Peridynamics.calc_force_density!(chunk::AbstractBodyChunk{<:MySystem}, t, Δt)`,
* `Peridynamics.log_system(::Type{<:MySystem}, options, dh)`.
