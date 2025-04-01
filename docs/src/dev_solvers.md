# Time solvers

!!! warning "Work in progress"
    This page is currently worked on.

## DynamicsRelaxation

## VelocityVerlet

## Custom time solver

define a type of subtype MySolver <: AbstractTimeSolver
define the init_time_solver!(vv::MySolver, dh::AbstractDataHandler) function
define the solve!(dh::AbstractDataHandler, vv::VelocityVerlet, options::AbstractJobOptions) function
