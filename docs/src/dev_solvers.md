# Time solvers

!!! warning "Work in progress"
    This page is currently worked on.

## VelocityVerlet

A time integration solver for the Velocity Verlet algorithm used for dynamic simulations.

## DynamicRelaxation

A time integration solver for the adaptive dynamic relaxation algorithm used for quasi-static simulations.

## Custom solvers

To create a custom time solver, the following steps are required:
- define a type `MySolver<:AbstractTimeSolver`. It has to be a subtype of the type `AbstractTimeSolver`.
- define the function `init_time_solver!(vv::MySolver, dh::AbstractDataHandler)`.
- define the function `solve!(dh::AbstractDataHandler, vv::VelocityVerlet, options::AbstractJobOptions)`.
