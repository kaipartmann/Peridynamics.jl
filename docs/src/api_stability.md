# API stability

Peridynamics.jl has three API tiers. The tier a name belongs to decides what you can rely
on, and every docstring says which tier it is in.

| Tier | How it is declared | How you write it | What it promises |
|:---|:---|:---|:---|
| **User API** | `export` | `Body`, `submit`, … | Stable. A change is breaking and is listed in [`NEWS.md`](https://github.com/kaipartmann/Peridynamics.jl/blob/main/NEWS.md). |
| **Extension API** | `public` | `Peridynamics.storage_type` | Stable within a minor release. A rename is breaking and is listed in `NEWS.md`. |
| **Internal** | neither | `Peridynamics.alloc_field` | Nothing. It can change in any release without notice. |

## The user API

Everything the package exports. This is what a simulation script is written with, and it is
listed on the [Public API](@ref) page. You get all of it with `using Peridynamics`.

## The extension API

The names you need to add something of your own: a material, a constitutive model, a damage
model or a set of point parameters. They are listed on the [Extension API](@ref) page.

These names are not exported, so they are written out in full:

```julia
using Peridynamics

function Peridynamics.force_density_point!(storage, system::Peridynamics.BondSystem,
                                           mat::MyMaterial, params, t, Δt, i)
    # ...
end
```

or imported explicitly, which reads better in a longer file:

```julia
using Peridynamics: BondSystem, each_bond_idx, get_vector_diff, update_add_vector!
```

Keeping them unexported is deliberate. It keeps `using Peridynamics` small in a simulation
script, and it makes every extension point visible as such when you read the code.

!!! note "While the package is at version 0.x"
    Julia's semantic versioning treats a minor bump of a `0.x` version as breaking, and this
    package uses that. So the extension API can still change in a `0.x` minor release. Every
    such change is written down in `NEWS.md` together with what to replace it with.

## What is deliberately not in the extension API

Some things are documented but internal on purpose, because they are expected to change:

- **The family-level stress hooks** of the correspondence materials, e.g.
  `calc_first_piola_kirchhoff!`, `rkc_stress_integral!` and `monomial`. Their argument
  lists differ per material family and are meant to be unified. Write a constitutive model
  against [`first_piola_kirchhoff`](@ref Peridynamics.first_piola_kirchhoff) instead, which
  is stable and works on every correspondence family.
- **Everything the macros generate for you**, e.g. `allowed_material_kwargs`,
  `req_storage_fields`, `point_param_type` and `log_param_property`. You never write these,
  so they are free to change.
- **The macro expansion machinery**, e.g. `storage_fields_expr` and
  `derive_storage_type_params`.
- **The parallelization layer**: the data handlers, the chunk handler, the body chunk and
  the halo exchange. This is what makes multithreading and MPI work without you doing
  anything, and it has to stay free to change for a GPU backend.
- **The time solver and system interfaces.** Writing a new time solver or a new system is
  possible, but those interfaces are not settled yet and will be stabilized in a later
  release.

If you find yourself needing one of these, please open an issue. Adding a name to the
extension API is not a breaking change, so it can be done in a patch release. Removing one
is, which is why the list starts narrow.

## Checking a name

From Julia 1.11 on, the tier of a name can be queried:

```julia-repl
julia> Base.isexported(Peridynamics, :submit)              # user API
true

julia> Base.ispublic(Peridynamics, :force_density_point!)  # extension API
true

julia> Base.ispublic(Peridynamics, :alloc_field)           # internal
false
```

The package is tested against a checked-in snapshot of both lists, so a name never enters or
leaves a tier by accident.
