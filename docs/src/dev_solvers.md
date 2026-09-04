# Time solvers

!!! warning "This interface is internal"
    Unlike materials, constitutive models and damage models, the time solver interface is
    not part of the [Extension API](@ref). The names below can change in any release. It is
    documented because it helps when reading the package, and it will be stabilized in a
    later release. If you are writing a solver of your own, please open an issue, so that
    the interface can be settled around a real use case.

## The solvers of this package

| solver | for | notes |
|:---|:---|:---|
| [`VelocityVerlet`](@ref) | dynamic simulations | explicit time integration, the time step is estimated from the point parameter `bc` |
| [`DynamicRelaxation`](@ref) | quasi-static simulations | adaptive dynamic relaxation. Its time step is a pseudo time step of one, so the total load applied by a velocity boundary condition over `n` steps is `v * n` |
| [`NewtonKrylov`](@ref) | static simulations | matrix-free Newton-Krylov. It evaluates the force density several times per time step, for the Jacobian-vector products and the line search, which is why a history-dependent constitutive model cannot be used with it |

## Which fields a solver needs

A time solver does not allocate storage fields itself. It says, per field, whether it needs
that field and at which extent, and leaves the number of rows and the element type to the
field shape declared in the storage. That is what
[`FullField`](@ref Peridynamics.FullField) and [`EmptyField`](@ref Peridynamics.EmptyField)
express:

```julia
init_field_solver(::VelocityVerlet, ::AbstractSystem, ::Val{:velocity}) = FullField()
init_field_solver(::AbstractTimeSolver, ::AbstractSystem, ::Val{:velocity}) = EmptyField()
```

A storage that declares `velocity::PointVector` therefore gets a full field under
`VelocityVerlet` and a zero-length array under every other solver, without the storage having
to know which solver it is used with. `FullField(extent)` overrules the declared extent,
which is how `NewtonKrylov` gets `b_int` with halo entries.

The fields a solver needs are declared once as a field block, which a storage then
`@inherit`s:

```julia
Peridynamics.@storage_fields VelocityVerletFields begin
    @lth position::PointVector{Float64}
    displacement::PointVector
    velocity::PointVector
    velocity_half::PointVector
    acceleration::PointVector
    b_int::PointVector
    b_ext::PointVector
end
```

The package ships [`VelocityVerletFields`](@ref Peridynamics.VelocityVerletFields),
[`DynamicRelaxationFields`](@ref Peridynamics.DynamicRelaxationFields) and
[`NewtonKrylovFields`](@ref Peridynamics.NewtonKrylovFields). A storage that should work
with all three inherits all three. The fields a given run does not need cost nothing.

## Contract

| what you define | signature | required | default or fallback |
|:---|:---|:---:|:---|
| the type | `struct MySolver <: Peridynamics.AbstractTimeSolver ... end` | yes | |
| the setup | `Peridynamics.init_time_solver!(solver::MySolver, dh::AbstractDataHandler)` | yes | |
| the time loop | `Peridynamics.solve!(dh::AbstractDataHandler, solver::MySolver, options)` | yes | |
| the fields the solver reads | `req_point_data_fields_timesolver(::Type{MySolver})`, `req_bond_data_fields_timesolver(::Type{MySolver})`, `req_data_fields_timesolver(::Type{MySolver})` | yes | this is the part of the storage contract that is checked when a [`Job`](@ref) is created, and a missing field throws a [`StorageContractError`](@ref Peridynamics.StorageContractError) |
| the fields only this solver needs | a block declared with [`@storage_fields`](@ref Peridynamics.@storage_fields), plus `init_field_solver` methods | yes, if there are any | |
| the simulation log | `Peridynamics.log_timesolver(options, solver::MySolver)` | yes | |
| the registration | `Peridynamics.register_solver!(MySolver)` | yes | see the caveat below |
| that the force density is evaluated more than once per step | [`supports_history_dependence(::MySolver)`](@ref Peridynamics.supports_history_dependence) `= false` | only then | `true`, see [Constitutive models](@ref) |

!!! warning "`register_solver!` is a package-wide commitment"
    It makes the new solver's field contract apply to every storage in the package,
    including the ones this package ships, so it is not something to do lightly.
    `NewtonKrylov` itself is deliberately left unregistered in
    `src/time_solvers/newton_krylov.jl` for this reason. The internal storages already
    support it, but user-defined ones may not.

See also [Storages](@ref), [Constitutive models](@ref), [Systems](@ref),
[Extension API](@ref).
