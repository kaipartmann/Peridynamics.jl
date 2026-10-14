# Point parameters

Point parameters are the material properties of a single point, i.e. what
[`material!`](@ref) assigns to a point set. They are declared with
[`@params`](@ref Peridynamics.@params), which generates the struct, the constructor that
reads the keywords of `material!`, the list of allowed keywords and the simulation log lines
from one list of declarations, so these cannot drift apart. The tutorial
[Writing your own material](@ref tutorial_custom_material) declares a set in full.

## Contract

| what you define | signature | required | default or fallback |
|:---|:---|:---:|:---|
| the point parameters of a material | `Peridynamics.@params MyMaterial struct MyPointParameters ... end` | yes, for a material | [`InterfaceError`](@ref Peridynamics.InterfaceError) |
| the same parameters as another material | `Peridynamics.@params MyMaterial BBPointParameters` | alternative to the first form | only a type declared with `@params` can be shared |
| a reusable block of declarations | `Peridynamics.@params_fields MyBlock begin ... end` | no | |
| the parameters of a constitutive model | `Peridynamics.@cm_params MyModel struct MyModelParameters ... end` | no | the model has no parameters |
| the parameters of a damage model | `Peridynamics.@dmg_params MyDamage struct MyDamageParameters ... end` | no | the model has no parameters |

## Skeleton

```julia
Peridynamics.@params MyMaterial struct MyPointParameters
    @inherit StandardParameters
    @log "initial yield stress" @kwarg sigma_y σy = Inf
    @log "hardening modulus" @kwarg hardening Hiso = 0.0
end
```

`material!(body; horizon, rho, E, nu, Gc, sigma_y=250.0, hardening=1000.0)` then works, the
two new keywords are accepted, everything else is rejected as a typo, and both appear in the
simulation log under the labels given.

## The declarations

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
| `cm_params::ConstitutiveParameters` | the place for the parameters of the constitutive model |
| `dmg_params::DamageParameters` | the place for the parameters of the damage model |

## One rule for every right-hand side

Take `@derived (; δb) = get_bond_horizon(δ; bond_horizon)`. It splits into three parts:

| part of the right-hand side | what may appear there | what it becomes |
|:---|:---|:---|
| before the `;`, here `δ` | the parameters declared above, and `mat`, and `model` inside a `@cm_params` or `@dmg_params` body | ordinary arguments of the call |
| after the `;`, here `bond_horizon` | names of `material!` keywords, in shorthand | the keywords `material!` accepts |
| the left-hand side, here `δb` | | the parameters produced |

That is the whole scoping rule, and it is why the order of the declarations matters. The
keywords written after the `;` are the ones `material!` accepts, so the allowed keywords
cannot disagree with the call that reads them. A keyword the user did not pass is not
forwarded, so the provider decides on its own whether it is required:

```julia
get_discretization_params(; horizon, rho)                    # both required
get_elastic_params(; E=nothing, nu=nothing, G=nothing, ...)  # any two of six
```

A body may declare any number of groups, in any order.

## Types

A parameter declared without a type follows the float type of the simulation. Inside a
definition `FT` stands for that type.

| declaration | type of the parameter |
|:---|:---|
| `rho` | the float type of the simulation |
| `n::Int` | `Int`, for every simulation |
| `C::SArray{NTuple{4,3},FT,4,81}` | follows the float type of the simulation through `FT` |

The generated struct is parametric in `FT` when any parameter follows it, and it has one
more type parameter per model marker field:

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
simulation. Point parameters are `isbits`, which is what lets them be captured by value in a
kernel, so an array-valued parameter is an `SArray` and never an `Array`.

## Reusing declarations with `@inherit`

- `@inherit` includes all declarations of a parameter block or of the point parameters of
  another material, spliced in at the position of the `@inherit`.
- Two `@inherit`s may contribute the same parameter only if they declare it identically.
- A declaration in the body overrides an inherited one in place, keeping its position. This
  is how the dual-horizon material halves the bond constant of the bond-based one:

  ```julia
  Peridynamics.@params DHBBMaterial struct DHBBPointParameters
      @inherit BBPointParameters
      @derived bc = 0.5 * 18 * K / (π * δ^4)
  end
  ```

- Everything the package ships that can be inherited is listed in
  [Blocks you can inherit](@ref), and the same table is one call away:

  ```julia
  Peridynamics.block_table(Peridynamics.DiscretizationParameters)  # what a block exposes
  Peridynamics.block_table(BBMaterial())                           # what a material accepts
  ```

- Names inside a definition are resolved in the module of the definition first and in
  `Peridynamics` second, so the blocks and point parameters of this package can be written
  unqualified.

Blocks of your own are declared with
[`@params_fields`](@ref Peridynamics.@params_fields):

```julia
Peridynamics.@params_fields PlasticityParameters begin
    @log "initial yield stress" @kwarg sigma_y σy = Inf
    @log "hardening modulus" @kwarg hardening Hiso = 0.0
end
```

## Sharing the parameters of another material

If the parameters of a material are exactly those of another one, the second form of
`@params` makes it use the same type:

```julia
Peridynamics.@params MyMaterial BBPointParameters
```

This is how [`GBBMaterial`](@ref) uses the parameters of [`BBMaterial`](@ref) and
[`CRMaterial`](@ref) those of [`CMaterial`](@ref). Only point parameters defined with
`@params` can be shared this way.

## The model markers

A constitutive model and a damage model bring parameters of their own, declared with
[`@cm_params`](@ref Peridynamics.@cm_params) and
[`@dmg_params`](@ref Peridynamics.@dmg_params) in the same language. The point parameters of
the material give them a place with a marker field:

| marker field | who fills it | resolved to |
|:---|:---|:---|
| `cm_params::ConstitutiveParameters` | the constitutive model of the material | the type declared with `@cm_params`, or `Nothing` |
| `dmg_params::DamageParameters` | the damage model of the material | the type declared with `@dmg_params`, or `Nothing` |

A material with a marker supports every model, parameterized or not, without knowing any of
them. The keywords a model declares are accepted by `material!` exactly when the body's
model reads them, and the model's parameters are read flat off the point parameters, e.g.
`params.Gc`. A declaration of a model may read every parameter declared above the marker,
which is why the marker comes last.

The standard fracture parameters `Gc` and `εc` are declared this way. They belong to
[`CriticalStretch`](@ref) and not to the material, and `@inherit StandardParameters`
includes the `dmg_params` marker for them.

See also [Materials](@ref), [Storages](@ref), [Damage models](@ref),
[Constitutive models](@ref), [Blocks you can inherit](@ref), [Extension API](@ref).
