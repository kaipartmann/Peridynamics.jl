# Developing custom time solvers

The following guide focuses on how developers can create custom time solvers with the package `Peridynamics.jl`.
In the future, support for [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/) may be added.

!!! warning "Private API"
    The following guide is a work in progress and may change in the future.
    All features documented here are part of the private API and may change without notice.

Since the private API is used, for readability, we import a few types and functions into the current scope:
```julia
using Peridynamics
using Peridynamics: AbstractTimeSolver, AbstractDataHandler
using Peridynamics: export_reference_results
using Base.Threads
```

## The time solver type

First, a custom time solver must be defined as a `struct` that is a subtype of the abstract type `AbstractTimeSolver`, e.g.:
```julia
struct MyCustomSolver <: AbstractTimeSolver
    n_steps::Int
    # other fields
end
```
The custom solver can have fields like a time step size, a time step counter, etc. that are needed for the solver.
The solver type can also be a `mutable struct` and contain mutable fields.

## Required time solver interface functions

The following interface functions *must* be implemented for the custom solver type:
  - `Peridynamics.init_time_solver!`
  - `Peridynamics.solve!`

The function `Peridynamics.init_time_solver!` is called to initialize the custom solver with the data handler and options:
```julia
function Peridynamics.init_time_solver!(solver::MyCustomSolver, dh::AbstractDataHandler, options)
    # initialize or change fields in `solver`
    # keep empty and just return nothing if the solver does not need initialization
    return nothing
end
```

Inside the function `Peridynamics.solve!`, all logic of the custom time solver is implemented.
When the `solve!` function is called, the data handler `dh` is fully initialized with the initial state and the custom solver was initialized with `Peridynamics.init_time_solver!`.
The data handler `dh` has to be returned.
For example, the function could look like this:
```julia
function Peridynamics.solve!(dh::AbstractDataHandler, solver::MyCustomSolver, options)
    # run the time solver and change data in `dh`
    export_reference_results(dh, options)
    for n in 1:solver.n_steps
        custom_timestep!(dh, options, n)
    end
    return dh
end
```
To export reference results, the function `export_reference_results` is called with the data handler and options.
The function `custom_timestep!` in this example should be defined for the custom solver type and should contain the logic for one time step.

!!! compat "Different data handlers for parallel simulations"
    The data handler `dh` should not be accessed directly in the `solve!` function, if the type is just `AbstractDataHandler`.
    Instead, custom methods like `custom_timestep!` should be defined using the types:
      - `Peridynamics.AbstractThreadsBodyDataHandler`: Multithreading simulations with only one body
      - `Peridynamics.AbstractThreadsMultibodyDataHandler`: Multithreading simulations with multiple bodies and contact
      - `Peridynamics.AbstractMPIBodyDataHandler`: MPI simulations with only one body
    The solver logic will be different regarding the data handler type, therefore custom methods should be defined for each data handler type.

!!! compat "Generic functions using chunks"
    The `Peridynamics.BodyChunk` type is the generic type for a chunk of data in all data handlers.
    If you write functions utilizing this type, you can reuse a lot of code and use multithreading or MPI parallelism with the same function.

To create a custom time step function for multithreading single-body simulations, the function could look like this:
```julia
function custom_timestep!(dh::AbstractThreadsBodyDataHandler, options, n)
    @threads :static for chunk in dh.chunks
        # each thread accesses a `chunk`
        apply_timestep!(chunk, n)
    end
    return nothing
end
```
The same functionality can be implemented for MPI simulations:
```julia
function custom_timestep!(dh::AbstractMPIBodyDataHandler, options, n)
    # each process has only one chunk
    apply_timestep!(dh.chunk, n)
    return nothing
end
```

## Optional time solver interface functions

The following interface functions *can* be implemented for the custom solver type, if needed:
  - `Peridynamics.init_field_solver`
  - `Peridynamics.log_timesolver`
  - `Peridynamics.req_point_data_fields_timesolver!`
  - `Peridynamics.req_bond_data_fields_timesolver!`
  - `Peridynamics.req_data_fields_timesolver!`
