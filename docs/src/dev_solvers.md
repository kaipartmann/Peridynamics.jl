# Time solvers

!!! warning "This interface is internal"
    Unlike materials, constitutive models and damage models, the time solver interface is
    not part of the [Extension API](@ref). The names below can change in any release. It is
    documented because it helps when reading the package, and it will be stabilized in a
    later release. If you are writing a solver of your own, please open an issue, so that
    the interface can be settled around a real use case.

## The solvers of this package

- [`VelocityVerlet`](@ref): explicit time integration for dynamic simulations.
- [`DynamicRelaxation`](@ref): adaptive dynamic relaxation for quasi-static simulations.
  Its time step is a pseudo time step of one, so the total load applied by a velocity
  boundary condition over `n` steps is `v * n`.
- [`NewtonKrylov`](@ref): a matrix-free Newton-Krylov solver for static simulations. It
  evaluates the force density several times per time step, for the Jacobian-vector
  products and for the line search, which is why a history-dependent constitutive model
  cannot be used with it.

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
`VelocityVerlet` and a zero-length array under every other solver, without the storage
having to know which solver it is used with. `FullField(extent)` overrules the declared
extent, which is how `NewtonKrylov` gets `b_int` with halo entries.

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

Independently of that, a solver declares the fields it reads, which is the part of the
storage contract that is checked when a [`Job`](@ref) is created:

- `req_point_data_fields_timesolver(::Type{MySolver})`
- `req_bond_data_fields_timesolver(::Type{MySolver})`
- `req_data_fields_timesolver(::Type{MySolver})`

A missing field then results in a `Peridynamics.StorageContractError` that names the field
and why it is required, instead of an error deep inside the time loop.

## Custom solvers

To add a time solver:

- define `MySolver <: Peridynamics.AbstractTimeSolver`,
- define `Peridynamics.init_time_solver!(solver::MySolver, dh::AbstractDataHandler)`,
- define `Peridynamics.solve!(dh::AbstractDataHandler, solver::MySolver, options)`,
- define the three `req_*_fields_timesolver` methods above,
- declare a field block with `@storage_fields` and add `init_field_solver` methods for the
  fields that only this solver needs,
- define `Peridynamics.log_timesolver(options, solver::MySolver)` for the simulation log,
- call `Peridynamics.register_solver!(MySolver)`.

If the solver evaluates the force density more than once per time step, also define

```julia
Peridynamics.supports_history_dependence(::MySolver) = false
```

so that a history-dependent constitutive model is rejected when the `Job` is created rather
than integrating its history several times per step.
