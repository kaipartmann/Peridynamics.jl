# Storages

A storage holds every field that changes during a simulation, one array per quantity. It is
declared with [`@storage`](@ref Peridynamics.@storage), and the same declaration language
also declares the state of a constitutive model and the state of a damage model. The
tutorial [Writing your own material](@ref tutorial_custom_material) declares a storage in
full.

## Contract

| what you define | signature | required | default or fallback |
|:---|:---|:---:|:---|
| the storage of a material | `Peridynamics.@storage MyMaterial struct MyStorage ... end` | yes, for a material | [`InterfaceError`](@ref Peridynamics.InterfaceError) |
| a reusable block of fields | `Peridynamics.@storage_fields MyBlock begin ... end` | no | |
| the state of a constitutive model | `Peridynamics.@cm_storage MyModel struct MyState ... end` | no | the model has no state |
| the state of a damage model | `Peridynamics.@dmg_storage MyDamage struct MyState ... end` | no | the model has no state |
| a field no shape describes | [`init_field(mat, solver, system, ::Val{:field})`](@ref Peridynamics.init_field) | for a container-typed field | a shaped field is allocated by its shape |

## Skeleton

```julia
Peridynamics.@storage MyMaterial struct MyStorage
    @inherit VelocityVerletFields DynamicRelaxationFields NewtonKrylovFields
    @inherit BondLengthCache
    my_point_field::PointScalar
    @htl my_accumulated_field::PointVector
    my_bond_field::BondSymTensor
    cm_state::ConstitutiveState
    dmg_state::DamageState
end
```

## Field shapes

A field shape says what a field *means*. It determines the container type, the element type
and how the field is allocated, so a shaped field needs no `init_field` method.

| shape | container | rows | entries |
|:---|:---|:---|:---|
| [`PointScalar`](@ref Peridynamics.PointScalar) | `Vector` | | points |
| [`PointVector`](@ref Peridynamics.PointVector) | `Matrix` | `get_n_dim(system)` | points |
| [`PointTensor`](@ref Peridynamics.PointTensor) | `Matrix` | `get_n_dim(system)^2` | points |
| [`PointSymTensor`](@ref Peridynamics.PointSymTensor) | `Matrix` | `d * (d + 1) ÷ 2` | points |
| [`PointField{N}`](@ref Peridynamics.PointField) | `Matrix` | `N` | points |
| [`BondScalar`](@ref Peridynamics.BondScalar) | `Vector` | | bonds |
| [`BondVector`](@ref Peridynamics.BondVector) | `Matrix` | `get_n_dim(system)` | bonds |
| [`BondTensor`](@ref Peridynamics.BondTensor) | `Matrix` | `get_n_dim(system)^2` | bonds |
| [`BondSymTensor`](@ref Peridynamics.BondSymTensor) | `Matrix` | `d * (d + 1) ÷ 2` | bonds |
| [`BondField{N}`](@ref Peridynamics.BondField) | `Matrix` | `N` | bonds |
| [`DofVector`](@ref Peridynamics.DofVector) | `Vector` | | degrees of freedom |

**A point shape makes the field point data**, which is what can be exported to VTK files. A
bond or dof shape is not point data.

Every shape takes an optional element type, e.g. `PointScalar{Bool}` or
`PointVector{Float64}`. Without one the field follows the float type of the simulation, with
one it keeps that element type for every simulation. This is why `position` is declared
`PointVector{Float64}`. Bond vectors are position differences and a smaller float type loses
them over a large domain.

## Halo annotations

Only the halo exchange is annotated, and both annotations give the field one entry per local
**and** halo point.

| annotation | entries | who is updated |
|:---|:---|:---|
| none | `get_n_loc_points(system)` | nothing is exchanged |
| [`@lth`](@ref Peridynamics.@lth) | `get_n_points(system)` | the halo entries are written from the chunk that owns the points, which is what `position` needs |
| [`@htl`](@ref Peridynamics.@htl) | `get_n_points(system)` | the halo entries are added into the local entries of the owning chunk, which is what an accumulated force density needs |

A material never has to do anything for this beyond the annotation. The exchange is the same
code under multithreading and under MPI.

## Initial value

| written | result |
|:---|:---|
| omitted | zeros |
| a `Number`, e.g. `= 1.0` or `= true` | every entry that value |
| a `LinearAlgebra.UniformScaling`, e.g. `= I` or `= 2I` | every column that tensor |

A `UniformScaling` requires a shape whose rows are the entries of a square tensor, e.g.
`PointTensor`.

## Anything a shape does not describe

- A field declared with a plain container type, e.g. `Vector{Float64}`, still needs an
  [`init_field`](@ref Peridynamics.init_field) method that allocates it. It is not point
  data.
- An `init_field` method is also the escape hatch for a **shaped** field. It is more
  specific than the generic fallback and therefore wins, so a field keeps its shape, and
  with it its type, its size and its export status, while being filled by hand.
- A field that a time solver works with, e.g. `velocity` or `residual`, has to be shaped. A
  solver says only whether it needs the field and at which extent and leaves the number of
  rows and the element type to the shape, so it cannot size a container-typed field.

## What you may read and write

A shaped field is one plain matrix with the quantity of a point or a bond in its columns, so
a kernel reads and writes whole columns as static vectors and tensors.

| function | returns | notes |
|:---|:---|:---|
| [`dims(x)`](@ref Peridynamics.dims) | the number of spatial dimensions as a `Val{N}` | `dims(system)`, `dims(storage)` or `dims(state)` |
| [`get_vector(M, i, dims)`](@ref Peridynamics.get_vector) | column `i` as an `SVector{N}` | |
| [`get_vector_diff(M, i, j, dims)`](@ref Peridynamics.get_vector_diff) | column `j` minus column `i` | the bond vector |
| [`update_vector!(M, i, v, dims)`](@ref Peridynamics.update_vector!) | writes, overwrites column `i` | |
| [`update_add_vector!(M, i, v, dims)`](@ref Peridynamics.update_add_vector!) | writes, adds to column `i` | how a force density accumulates |
| [`get_tensor(M, i, dims)`](@ref Peridynamics.get_tensor) | column `i` as an `SMatrix{N,N}` | column-major order |
| [`update_tensor!(M, i, A, dims)`](@ref Peridynamics.update_tensor!) | writes an `SMatrix{N,N}` into column `i` | |
| [`get_sym_tensor(M, i, dims)`](@ref Peridynamics.get_sym_tensor) | column `i` of a symmetric shape as an `SMatrix{N,N}` | Voigt order in the array |
| [`update_sym_tensor!(M, i, A, dims)`](@ref Peridynamics.update_sym_tensor!) | writes a symmetric `SMatrix{N,N}` into column `i` | |

The `Val{N}` is a compile time constant, so the call folds into plain indexing. A value
written back has to be a static vector or tensor of that same `N`, anything else is a
`MethodError`.

```julia
Δxij = get_vector_diff(storage.position, i, j, dims(system))
update_add_vector!(storage.b_int, i, b, dims(system))
εp = get_sym_tensor(state.bond_plastic_strain, idx, dims(storage))
```

## The generated type

The generated struct carries the number of spatial dimensions of the simulation, is generic
in its float type and is parametric in the array type of every field:

```julia
struct BBStorage{N,FT<:Real,M_F64<:AbstractMatrix{Float64},M_FT<:AbstractMatrix{FT},
                 V_FT<:AbstractVector{FT},DMS} <: AbstractStorage
    position::M_F64
    displacement::M_FT
    ...
    dmg_state::DMS
end
```

The macro derives these parameters from the field declarations, one per distinct combination
of element type and number of array dimensions, plus one unbounded parameter per nested
state (`cm_state`, `dmg_state`) that the storage declares, so a storage must not declare
type parameters of its own.

- `Peridynamics.storage_type(mat)` returns the instantiation with the arrays of the CPU,
  `storage_type(mat, Float32)` the one with `Float32` arrays, and
  `storage_type(mat, Float64, Val(2))` the two-dimensional one.
- `Adapt.adapt(backend, storage)` moves a whole storage to another array backend.
- `N` comes first and is always there, so [`get_n_dim`](@ref Peridynamics.get_n_dim) and
  with it [`dims`](@ref Peridynamics.dims) work on a storage exactly as they do on a system.
  The storage of a body chunk is built with the `N` of its system, so the two can never
  disagree, which is what lets a constitutive model hook name its dimension without a system
  in scope.
- Dispatch on the storage *type* therefore has to be written `::Type{<:MyStorage}` instead
  of `::Type{MyStorage}`, while dispatch on a storage *value*, e.g. `::MyStorage`, is
  unchanged.

## Reusing fields with `@inherit`

- `@inherit` includes all fields of another storage or of a field block, spliced in at the
  position of the `@inherit`, keeping their order.
- Several `@inherit`s may contribute the same field only if they declare it identically.
- A field declared in the body overrides an inherited field of the same name and keeps its
  position, which is how a family changes an annotation, e.g. from `b_int::PointVector` to
  `@htl b_int::PointVector`.
- Every field block and every storage this package ships is listed in
  [Blocks you can inherit](@ref), and `Peridynamics.block_table(VelocityVerletFields)`
  prints the same table.
- Blocks of your own are declared with
  [`@storage_fields`](@ref Peridynamics.@storage_fields), by a top-level statement earlier
  than the storage that inherits them:

  ```julia
  Peridynamics.@storage_fields MyFamilyFields begin
      my_point_field::PointScalar
      my_bond_field::BondScalar
  end
  ```

## `BondLengthCache`

- [`BondLengthCache`](@ref Peridynamics.BondLengthCache) is an optional block carrying
  `bond_length`, the current length of every bond, filled once per point and per time step
  before the damage model and the force density run.
- It is worth 8 bytes per bond for a material whose force density needs the current length
  anyway, e.g. `BBMaterial` and `OSBMaterial`, and it is not worth it for one that does not,
  e.g. `CMaterial`.
- Either way [`current_bond_length`](@ref Peridynamics.current_bond_length) and
  [`bond_stretch`](@ref Peridynamics.bond_stretch) read the length, and a kernel does not
  have to know which of the two it is. Nothing reaches `storage.bond_length` directly, not
  even the materials of this package.

## The state markers

| marker field | who fills it | resolved to |
|:---|:---|:---|
| `cm_state::ConstitutiveState` | the constitutive model of the material | the type declared with `@cm_storage`, or `Nothing` |
| `dmg_state::DamageState` | the damage model of the material | the type declared with `@dmg_storage`, or `Nothing` |

Each marker contributes one type parameter, which
[`storage_type`](@ref Peridynamics.storage_type) fills with the state of the model that is
actually used, so the storage stays concrete and a model without state costs a zero-size
field. The fracture bookkeeping is therefore not declared by the material. It belongs to the
damage model, which carries [`BondFracFields`](@ref Peridynamics.BondFracFields) in its
state. A material never touches those fields by name, it asks
[`bond_is_active`](@ref Peridynamics.bond_is_active) and
[`get_damage`](@ref Peridynamics.get_damage) instead.

## The storage contract

A storage has to contain every field that is read by the code it inherits, and which fields
these are is not all known when `@storage` is expanded. The complete contract is checked
once when a [`Job`](@ref) is created.

| who requires fields | example | error |
|:---|:---|:---|
| the material family | every material of the RKC family needs the fields of [`RKCFields`](@ref Peridynamics.RKCFields) | [`StorageContractError`](@ref Peridynamics.StorageContractError) |
| the damage model | a model with a state needs the field `dmg_state` | [`StorageContractError`](@ref Peridynamics.StorageContractError) |
| the constitutive model | a history-dependent model needs the field `cm_state` | [`HistoryDependenceError`](@ref Peridynamics.HistoryDependenceError) |
| the time solver | [`NewtonKrylov`](@ref) needs `residual`, `Δu` and further buffers | [`StorageContractError`](@ref Peridynamics.StorageContractError) |

The fracture bookkeeping is not part of the contract. Everything outside the damage model
reads it through `bond_is_active` and `get_damage`, which behave neutrally for a model that
carries none.

See also [Materials](@ref), [Point parameters](@ref), [Damage models](@ref),
[Constitutive models](@ref), [Time solvers](@ref), [Blocks you can inherit](@ref),
[Extension API](@ref).
