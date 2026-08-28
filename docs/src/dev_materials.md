# Materials

A material decides which peridynamic formulation a body is simulated with. This page is the
manual of the declaration language a material, a damage model and a constitutive model are
written in. The complete, runnable examples are the three tutorials
[Writing your own material](@ref tutorial_custom_material),
[Writing your own damage model](@ref tutorial_custom_damage_model) and
[Writing your own constitutive model](@ref tutorial_custom_constitutive_model), and
[Extension API](@ref) lists every name used here.

## The materials of this package

- [`BBMaterial`](@ref): bond-based peridynamics.
- [`DHBBMaterial`](@ref): dual-horizon bond-based peridynamics.
- [`GBBMaterial`](@ref): generalized bond-based peridynamics.
- [`OSBMaterial`](@ref): ordinary state-based peridynamics, also called linear peridynamic
  solid (LPS).
- [`CMaterial`](@ref): the correspondence formulation.
- [`CRMaterial`](@ref): the correspondence formulation with stress rotation.
- [`RKCMaterial`](@ref): reproducing kernel peridynamics with bond-associated higher order
  integration.
- [`RKCRMaterial`](@ref): the reproducing kernel formulation with stress rotation.
- [`BACMaterial`](@ref): the bond-associated correspondence formulation of Chen and Spencer.
- [`CKIMaterial`](@ref): continuum-kinematics-inspired peridynamics.

## What a material consists of

A material needs four things:

1. **a type**, whose supertype says which system it is discretized on,
2. **its point parameters**, declared with [`@params`](@ref Peridynamics.@params),
3. **its storage**, declared with [`@storage`](@ref Peridynamics.@storage),
4. **the force density calculation**, a method of
   [`force_density_point!`](@ref Peridynamics.force_density_point!).

Everything these need is part of the [Extension API](@ref), so it is written as
`Peridynamics.<name>` or imported explicitly. Inside the force density a material walks the
bonds of a point with [`each_bond_idx`](@ref Peridynamics.each_bond_idx), reads the bond off
`system.bonds`, and reads the storage through the fields of the blocks it inherited and its
own fields:

```julia
function Peridynamics.force_density_point!(storage::MyStorage, system::BondSystem,
                                           mat::MyMaterial, params, t, Δt, i)
    (; bonds, correction, volume) = system
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(storage.position, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        ω = storage.bond_active[bond_id] * surface_correction_factor(correction, bond_id)
        b = ω * params.bc * ε * volume[j] / l .* Δxij
        update_add_vector!(storage.b_int, i, b)
    end
    return nothing
end
```

Which bonds are broken was decided right before by the damage model, so the force density
multiplies `bond_active` in and never changes it.

### Two names that are not free

A material on a bond system is reached by name in two places, so these two names have to be
spelled exactly like this:

- **The material needs a field `dmgmodel`**, because every bond system material is asked
  for its damage model before the force density is evaluated. Give it a type parameter and
  a keyword, as every material of this package does:
  ```julia
  struct MyMaterial{Correction,DM} <: Peridynamics.AbstractBondSystemMaterial{Correction}
      dmgmodel::DM
  end
  ```
- **The point parameters need a `bc`**, the bond constant, because the stable time step of
  an explicit solver is estimated from it. `@inherit StandardParameters` already derives it.
  A material that derives its own has to keep the name. If the bond stiffness is not
  constant over the family, declare `bc` as its largest value, so that the estimate stays
  on the safe side. Such a material also defines
  [`critical_stretch`](@ref Peridynamics.critical_stretch) and
  [`energy_release_rate`](@ref Peridynamics.energy_release_rate), because the relation
  between `Gc` and `εc` depends on the micro-modulus.

## Declaring the point parameters

Point parameters are the material properties of a single point: what [`material!`](@ref)
assigns to a point set. [`@params`](@ref Peridynamics.@params) generates the struct, the
constructor that reads the keywords of `material!`, the list of allowed keywords and the
simulation log lines from one list of declarations, so they cannot drift apart:

```julia
Peridynamics.@params MyMaterial struct MyPointParameters
    @inherit StandardParameters
    @log "initial yield stress" @kwarg sigma_y σy = Inf
    @log "hardening modulus" @kwarg hardening Hiso = 0.0
end
```

`material!(body; horizon, rho, E, nu, Gc, sigma_y=250.0, hardening=1000.0)` then works,
the two new keywords are accepted, everything else is rejected as a typo, and both appear
in the simulation log under the labels given.

### The declarations

| declaration | meaning |
|:---|:---|
| `rho` | required keyword `rho` |
| `sigma_y = Inf` | keyword `sigma_y`, defaulting to `Inf` |
| `C1 = 30 * μ / (π * δ^4)` | keyword `C1`, defaulting to an expression of the parameters above |
| `n::Int = 4` | keyword pinned to a concrete type |
| `@kwarg gamma_c gammac = 1e-10` | keyword `gamma_c`, parameter `gammac` |
| `@derived bc = 18 * K / (π * δ^4)` | computed, **not** a keyword |
| `@derived (; δ, rho) = get_discretization_params(; horizon, rho)` | a group of parameters computed by one call |
| `@log "shear modulus" G` | also write `G` to the simulation log |
| `@log "yield stress" sigma_y = Inf` | declare and log in one line |
| `@inherit StandardParameters` | include the declarations of another block |
| `dmg_params::DamageParameters` | the place for the parameters of the damage model |

### One rule for every right-hand side

```
@derived (; δb) = get_bond_horizon(δ; bond_horizon)
            ~~                     ~  ~~~~~~~~~~~~
      parameters produced   before `;`: the         after `;`: names of
                            parameters declared     `material!` keywords,
                            above, and `mat`        in shorthand
```

That is the whole scoping rule, and it is why the order of the declarations matters. The
keywords written after `;` are the ones `material!` accepts, so the allowed keywords cannot
disagree with the call that reads them. A keyword the user did not pass is not forwarded,
so the provider decides on its own whether it is required:

```julia
get_discretization_params(; horizon, rho)                    # both required
get_elastic_params(; E=nothing, nu=nothing, G=nothing, ...)  # any two of six
```

A body may declare any number of groups, in any order.

### Types

A parameter declared without a type follows the float type of the simulation. Inside a
definition `FT` stands for that type, so a parameter whose type is built from it stays
generic as well:

| declaration | type of the parameter |
|:---|:---|
| `rho` | the float type of the simulation |
| `n::Int` | `Int`, for every simulation |
| `C::SArray{NTuple{4,3},FT,4,81}` | follows the float type of the simulation through `FT` |

The generated struct is parametric in `FT` when any parameter follows it, and it has one
more type parameter per marker field of a model (see below):

```julia
struct BBPointParameters{FT<:Real,DMP} <: AbstractPointParameters
    δ::FT
    rho::FT
    ...
    bc::FT
    dmg_params::DMP
end
```

`Peridynamics.point_param_type(mat)` returns the instantiation for the float type of the
simulation. Point parameters are `isbits`, which is what lets them be captured by value in
a kernel, so an array-valued parameter is an `SArray`, never an `Array`.

### Reusing parameters with `@inherit`

`@inherit` includes all declarations of a parameter block or of the point parameters of
another material. Two `@inherit`s may contribute the same parameter only if they declare it
identically, and a declaration in the body overrides an inherited one in place. This is how
the dual-horizon material halves the bond constant of the bond-based one:

```julia
Peridynamics.@params DHBBMaterial struct DHBBPointParameters
    @inherit BBPointParameters
    @derived bc = 0.5 * 18 * K / (π * δ^4)
end
```

Everything this package ships that can be inherited is listed in
[Blocks you can inherit](@ref), and the reference entry of every block says what it
exposes. The same table is one call away:

```julia
Peridynamics.block_table(Peridynamics.DiscretizationParameters)  # what a block exposes
Peridynamics.block_table(BBMaterial())                           # what a material accepts
```

Own blocks are declared with [`@params_fields`](@ref Peridynamics.@params_fields):

```julia
Peridynamics.@params_fields PlasticityParameters begin
    @log "initial yield stress" @kwarg sigma_y σy = Inf
    @log "hardening modulus" @kwarg hardening Hiso = 0.0
end
```

Names inside a definition are resolved in the module of the definition first and in
`Peridynamics` second, so the blocks and point parameters of this package can be written
unqualified.

### Sharing the parameters of another material

If the parameters of a material are exactly those of another one, the second form of
`@params` makes it use the same type:

```julia
Peridynamics.@params MyMaterial BBPointParameters
```

This is how [`GBBMaterial`](@ref) uses the parameters of [`BBMaterial`](@ref) and
[`CRMaterial`](@ref) those of [`CMaterial`](@ref). Only point parameters defined with
`@params` can be shared this way.

### Parameters of the models

A constitutive model or a damage model can bring parameters of its own, declared with
[`@cm_params`](@ref Peridynamics.@cm_params) or
[`@dmg_params`](@ref Peridynamics.@dmg_params) in the same language. The point parameters of
the material give them a place with a marker field:

```julia
Peridynamics.@params MyMaterial struct MyPointParameters
    @inherit DiscretizationParameters ElasticParameters
    @derived bc = 18 * K / (π * δ^4)
    cm_params::ConstitutiveParameters
    dmg_params::DamageParameters
end
```

The type behind a marker is resolved per model when the point parameter type is
instantiated, and it is `Nothing` for a model without parameters, so a material with the
marker supports every model, parameterized or not, without knowing any of them. The
keywords a model declares are accepted by `material!` exactly when the body's model reads
them, and the model's parameters are read flat off the point parameters, e.g. `params.Gc`.
A declaration of a model may read every parameter declared above the marker, which is why
the marker comes last.

The standard fracture parameters `Gc` and `εc` are declared this way: they belong to
[`CriticalStretch`](@ref), not to the material, and `@inherit StandardParameters` includes
the `dmg_params` marker for them.

## Declaring the storage

A storage holds every field that changes during a simulation, one array per quantity. It is
declared with [`@storage`](@ref Peridynamics.@storage). Every field is declared either with
a **field shape** or with a plain container type, and it can carry one of the halo
annotations `@lth` or `@htl`.

A field shape says what a field *means*. It determines the container type, the element type
and how the field is allocated, so a shaped field needs no `init_field` method:

| shape | container | rows | entries |
|:---|:---|:---|:---|
| `PointScalar` | `Vector` | – | points |
| `PointVector` | `Matrix` | `get_n_dim(system)` | points |
| `PointTensor` | `Matrix` | `get_n_dim(system)^2` | points |
| `PointSymTensor` | `Matrix` | `d * (d + 1) ÷ 2` | points |
| `PointField{N}` | `Matrix` | `N` | points |
| `BondScalar` | `Vector` | – | bonds |
| `BondVector` | `Matrix` | `get_n_dim(system)` | bonds |
| `BondTensor` | `Matrix` | `get_n_dim(system)^2` | bonds |
| `BondSymTensor` | `Matrix` | `d * (d + 1) ÷ 2` | bonds |
| `BondField{N}` | `Matrix` | `N` | bonds |
| `DofVector` | `Vector` | – | degrees of freedom |

**A point shape makes the field point data**, which is what can be exported to VTK files. A
bond or dof shape is not point data.

Every shape takes an optional element type, e.g. `PointScalar{Bool}` or
`PointVector{Float64}`. Without one the field follows the float type of the simulation,
with one it keeps that element type for every simulation. This is why `position` is
declared `PointVector{Float64}`: bond vectors are position differences and a smaller float
type loses them over a large domain.

Only the halo exchange is annotated, and both annotations give the field one entry per
local *and* halo point:

- `@lth`: local-to-halo exchange, so the halo entries are updated from the chunk that owns
  the points.
- `@htl`: halo-to-local exchange, so the halo entries are added to the local entries of
  the owning chunk.

The initial value is given with `= value`, where `value` is a `Number` that fills every
entry or a `LinearAlgebra.UniformScaling` such as `I` or `2I` that writes that tensor into
every column of a square tensor shape. It defaults to zero.

A field declared with a container type, e.g. `Vector{Float64}`, still needs an
[`init_field`](@ref Peridynamics.init_field) method that allocates it. This is the right
choice for anything a shape cannot express. The one exception is a field that a time solver
works with, e.g. `velocity` or `residual`. A solver says only whether it needs the field and
at which extent and leaves the number of rows and the element type to the shape, so such a
field has to be shaped.

An `init_field` method is also the escape hatch for a shaped field. It is more specific
than the generic fallback and therefore wins, so a field can keep its shape, and with it
its type, its size and its export status, while being filled by hand.

### The generated type

The generated struct is parametric in the array type of every field and generic in the
float type of the simulation:

```julia
struct BBStorage{FT<:Real,M_F64<:AbstractMatrix{Float64},M_FT<:AbstractMatrix{FT},
                 V_FT<:AbstractVector{FT},V_Int<:AbstractVector{Int},
                 V_Bool<:AbstractVector{Bool}} <: AbstractStorage
    position::M_F64
    displacement::M_FT
    ...
end
```

The macro derives these parameters from the field declarations, one per distinct
combination of element type and number of dimensions, so a storage must not declare type
parameters of its own. `Peridynamics.storage_type(mat)` returns the instantiation with the
arrays of the CPU, `Peridynamics.storage_type(mat, Float32)` the one with `Float32` arrays,
and `Adapt.adapt(backend, storage)` moves a whole storage to another array backend.

Dispatch on the storage *type* therefore has to be written `::Type{<:MyStorage}` instead of
`::Type{MyStorage}`, while dispatch on a storage *value*, e.g. `::MyStorage`, is unchanged.

### Reusing fields with `@inherit`

`@inherit` includes all fields of another storage or of a field block. This is how a custom
material reuses the fields of the family it builds on, instead of copying them:

```julia
Peridynamics.@storage MyMaterial struct MyStorage
    @inherit RKCStorage
    my_point_field::PointScalar
    my_bond_field::BondScalar
end
```

The inherited fields keep their order and are spliced in at the position of the `@inherit`.
Several `@inherit`s may contribute the same field as long as they declare it identically. A
field declared in the body itself overrides an inherited field of the same name and keeps
its position, which is how a family changes an annotation, e.g. from `b_int::PointVector`
to `@htl b_int::PointVector`.

Every field block and every storage this package ships is listed in
[Blocks you can inherit](@ref), its reference entry shows the fields it exposes, and
`Peridynamics.block_table(VelocityVerletFields)` prints the same table. A storage that supports all three time solvers inherits their three
blocks:

```julia
Peridynamics.@storage BBMaterial struct BBStorage
    @inherit VelocityVerletFields DynamicRelaxationFields NewtonKrylovFields BondFracFields
    strain_energy_density::PointScalar
    dmg_state::DamageState
end
```

The marker `dmg_state::DamageState` is the place for the state of the damage model, see
below. Every storage of this package declares it, so every damage model runs on every
material.

Own field blocks are defined with [`@storage_fields`](@ref Peridynamics.@storage_fields):

```julia
Peridynamics.@storage_fields MyFamilyFields begin
    my_point_field::PointScalar
    my_bond_field::BondScalar
end
```

A block has to be defined by an earlier top-level statement than the storage that inherits
it.

## The storage contract

A storage has to contain every field that is read by the code it inherits. Which fields
these are depends on three things that are not all known when `@storage` is expanded:

1. the material family, e.g. every material of the RKC family needs the fields of
   `RKCFields`,
2. the damage model, e.g. a model with a state of its own needs the field `dmg_state`,
3. the time solver, e.g. `NewtonKrylov` needs `residual`, `Δu` and further buffers.

Therefore `@storage` only checks the part of the contract that follows from the material
type. The complete contract is checked once when a [`Job`](@ref) is created, and a missing
field results in a
[`StorageContractError`](@ref Peridynamics.StorageContractError) that names the field and
the reason why it is required.

## Constitutive models

The correspondence families ([`CMaterial`](@ref), [`RKCMaterial`](@ref),
[`BACMaterial`](@ref)) do not fix the stress-strain relation. They take a **constitutive
model** and ask it for the first Piola-Kirchhoff stress that belongs to a deformation
gradient, so a new material behavior usually does not need a new material at all. A model
is a subtype of [`AbstractConstitutiveModel`](@ref Peridynamics.AbstractConstitutiveModel)
that defines [`first_piola_kirchhoff`](@ref Peridynamics.first_piola_kirchhoff):

```julia
struct MyModel <: Peridynamics.AbstractConstitutiveModel end

function Peridynamics.first_piola_kirchhoff(::MyModel, storage, params, F)
    return ...
end
```

`RKCMaterial(model=MyModel())` then works, and so does every other family, on threads and
with MPI. The tutorial [Writing your own constitutive model](@ref
tutorial_custom_constitutive_model) writes a hyperelastic and a plastic model in full.

### Parameters of a model

A model that needs parameters of its own declares them with
[`@cm_params`](@ref Peridynamics.@cm_params), in the same language as `@params`. They
become keywords of `material!` for every material whose point parameters carry the marker
`cm_params::ConstitutiveParameters`, which the correspondence families do, and they are
read flat off the point parameters, e.g. `params.sigma_y`. A declaration may read every
material parameter declared above the marker, e.g. the shear modulus `μ`, the model
instance is available as `model` and the material as `mat`.

### History-dependent models

A model that integrates an internal state over time, such as plasticity, viscoelasticity or
creep, declares that state with [`@cm_storage`](@ref Peridynamics.@cm_storage), which
accepts the same field declarations as [`@storage`](@ref Peridynamics.@storage):

```julia
Peridynamics.@cm_storage MyPlasticModel struct MyPlasticState
    bond_plastic_strain::BondSymTensor
    bond_eqps::BondScalar
end
```

The state is reached inside the stress update with
[`constitutive_state`](@ref Peridynamics.constitutive_state), and the stress update then
takes two more arguments, the index of the evaluated quantity and the time step:

```julia
function Peridynamics.first_piola_kirchhoff(::MyPlasticModel, storage, params, F, idx, Δt)
    state = Peridynamics.constitutive_state(storage)
    ...
end
```

A model that needs no state defines the four-argument form above, which is bridged to this
one. What `idx` indexes follows from the material family, so the state is declared with the
matching field shapes:

| material family | `idx` | state shapes |
|:---|:---|:---|
| [`CMaterial`](@ref) | point index | `Point...` |
| [`RKCMaterial`](@ref), [`BACMaterial`](@ref) | bond index | `Bond...` |

The state has to be carried by the storage of the material, which the storages of the three
families above already do with the declaration `cm_state::ConstitutiveState`. It
contributes one type parameter to the storage, which
[`storage_type`](@ref Peridynamics.storage_type) fills with the state of the model that is
actually used, so the storage stays concrete and a model without state costs a zero-size
field. The state is chunk-local and is never exchanged between chunks, which is why the halo
annotations are not allowed in a `@cm_storage` definition.

### What is checked

Declaring a state makes a model history dependent, see
[`is_history_dependent`](@ref Peridynamics.is_history_dependent). A history-dependent
model may only be run by a time solver that evaluates the force density **once** per time
step. [`NewtonKrylov`](@ref) evaluates it several times, for the Jacobian-vector products
and the line search, so the model would integrate its history several times per step. This
is checked once when a [`Job`](@ref) is created and results in a
[`HistoryDependenceError`](@ref Peridynamics.HistoryDependenceError) that names the reason.

### Energy

[`strain_energy_density`](@ref Peridynamics.strain_energy_density) also takes the index, but
not the time step, and it must not change the state. It is called when a field is exported,
that is outside of the time integration.

## Damage models

A damage model decides when a bond fails. It is a subtype of
[`AbstractDamageModel`](@ref Peridynamics.AbstractDamageModel), and the one method it has to
define is [`calc_failure!`](@ref Peridynamics.calc_failure!), which is called once per local
point and per time step, right before the force density, with the time and the time step:

```julia
function Peridynamics.calc_failure!(storage, system, mat, ::MyDamage, paramsetup, t, Δt, i)
    (; εc) = get_params(paramsetup, i)
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ...
        storage.n_active_bonds[i] += storage.bond_active[bond_id]
    end
    return nothing
end
```

A method deactivates the bonds that fail, counts the ones that are still active in
`n_active_bonds`, and never breaks a bond whose `fail_permit` is `false`, because that is
how [`no_failure!`](@ref) and the pre-cracks are honored. The tutorial
[Writing your own damage model](@ref tutorial_custom_damage_model) writes a model with a
delay in full.

### Fracture parameters

A damage model owns its fracture parameters and declares them with
[`@dmg_params`](@ref Peridynamics.@dmg_params). Inheriting the
[`FractureParameters`](@ref Peridynamics.FractureParameters) block brings the standard pair
`Gc` and `εc`, resolved from the keywords `Gc` and `epsilon_c` of `material!` by
[`get_frac_params`](@ref Peridynamics.get_frac_params), and a model adds keywords of its
own next to it:

```julia
Peridynamics.@dmg_params MyDamage struct MyDamageParameters
    @inherit FractureParameters
    @log "failure delay" @kwarg tau τ
end
```

Two things come for free with the block. The conversion between `Gc` and `εc` is the
default of every damage model, and it goes through
[`critical_stretch`](@ref Peridynamics.critical_stretch) and
[`energy_release_rate`](@ref Peridynamics.energy_release_rate), which dispatch on the damage
model **and** the material, because the relation follows from the micro-modulus. A material
with another micro-modulus defines those two, always both, and never `get_frac_params`. And
[`has_fracture`](@ref Peridynamics.has_fracture), which decides whether the bonds of a point
set may fail at all, reads `Gc` and `εc` by default, so leaving both keywords out switches
fracture off as it does for `CriticalStretch`. A model that reads *other* keywords defines
`get_frac_params`, a model whose own parameters are the fracture parameters defines
`has_fracture`.

A material carries the parameters of its damage model in the `dmg_params::DamageParameters`
marker field, which `@inherit StandardParameters` includes.

### State of its own

A model that needs per-bond variables, e.g. an accumulated damage, declares them with
[`@dmg_storage`](@ref Peridynamics.@dmg_storage), which accepts the same field declarations
as [`@storage`](@ref Peridynamics.@storage):

```julia
Peridynamics.@dmg_storage MyDamage struct MyDamageState
    bond_damage::BondScalar
end
```

The state is reached with [`damage_state`](@ref Peridynamics.damage_state), and a material
carries it by declaring the field `dmg_state::DamageState`, which every storage of this
package does, so every shipped material takes a damage model of yours. A material that declares that field supports **every** damage model, stateful
or not, without knowing any of them. The model brings its own arrays instead of the material
having to allocate them for it. A model without state answers `nothing`, and no arrays are
allocated at all.

Unlike a constitutive state, a damage state does not make anything history dependent. A
damage model advances its state in `calc_failure!`, which every time solver calls exactly
once per step, so a stateful damage model stays usable under [`NewtonKrylov`](@ref).

### Softening a bond instead of deleting it

Deleting a bond is a jump in the moment matrix of a reproducing kernel material, and no
regularization of its inverse can absorb a jump in its input. A model can therefore let a
bond fade out instead, through two hooks that both default to one:

| hook | scales |
|---|---|
| [`bond_integrity`](@ref Peridynamics.bond_integrity) | the stress and the strain energy the bond carries |
| [`kinematic_weight`](@ref Peridynamics.kinematic_weight) | what the bond contributes to the moment matrix and the gradient weights |

A material says with [`supports_bond_integrity`](@ref Peridynamics.supports_bond_integrity)
and [`supports_kinematic_weight`](@ref Peridynamics.supports_kinematic_weight) whether its
force path calls the hooks. [`RKCMaterial`](@ref) and [`RKCRMaterial`](@ref) do. Combining
a softening model with a material that ignores the hooks fails once when the `Job` is
created, instead of silently not softening.

A model that softens also defines [`calc_damage!`](@ref Peridynamics.calc_damage!), because
the default damage of a point is the fraction of deleted bonds, which is not what a
softening model means.

## Exporting fields of your own

Every point field of a storage can be named in the `fields` keyword of a [`Job`](@ref) and
is written to the VTK files as it is. A quantity that is not a storage field, or a bond
field that has to be reduced to one value per point, is exported through
[`export_field`](@ref Peridynamics.export_field), and its name is announced with
[`custom_field`](@ref Peridynamics.custom_field), so that asking for it is not rejected as a
typo:

```julia
Peridynamics.custom_field(::Type{<:MyStorage}, ::Val{:bond_damage_avg}) = true

function Peridynamics.export_field(::Val{:bond_damage_avg}, mat, system, storage::MyStorage,
                                   paramsetup, t)
    n = Peridynamics.get_n_loc_points(system)
    out = zeros(n)
    for i in 1:n
        bond_ids = Peridynamics.each_bond_idx(system, i)
        out[i] = sum(@view storage.bond_damage[bond_ids]) / length(bond_ids)
    end
    return out
end
```

The returned array has one entry per local point.
