# [Manual](@id manual)

!!! tip "Typical workflow"
    `Peridynamics.jl` simulations follow this short workflow:

    1. Create a job
    2. Submit the job
   

## How to create a job

Currently two different simulation job types are available:

- [`PDSingleBodyAnalysis`](@ref)
  - Define the spatial discretization - create a [`PointCloud`](@ref)
  - Define the material, see [`BondBasedMaterial`](@ref)
  - Define the temporal discretization with [`TimeDiscretization`](@ref)
  - Define export settings with [`ExportSettings`](@ref)
  - Add predefined cracks with [`PreCrack`](@ref)
  - Add boundary conditions, e.g. [`VelocityBC`](@ref), [`ForceDensityBC`](@ref), [`PosDepVelBC`](@ref)
  - Add initial conditions, e.g. [`VelocityIC`](@ref)
- [`PDContactAnalysis`](@ref)
  - Create a [`BodySetup`](@ref) for every body (similar to `PDSingleBodyAnalysis`)
  - Define the [`Contact`](@ref) between the different bodies
  - Define the temporal discretization with [`TimeDiscretization`](@ref)
  - Define export settings with [`ExportSettings`](@ref)
  

## How to submit a job

Use the [`submit`](@ref) function:  
```julia
submit(job)
```

## Point clouds

There are three possible ways of creating [`PointCloud`](@ref)s.

1. Read `.inp` files with the [`read_inp`](@ref) function
2. Create a rectangular [`PointCloud`](@ref) by specifying the edge lengths and the point spacing
3. Assemble point clouds manually by specifying point coordinates and point volume


## Material models

Currently only the [`BondBasedMaterial`](@ref) is implemented.
More material models will be added in the future!


## Temporal discretization

The temporal discretization is controlled with the [`TimeDiscretization`](@ref) type.
By default, explicit time integration with the Velocity Verlet algorithm is used.
Quasistatic simulations can currently be performed using the adaptive dynamic relaxation algorithm by [Kilic and Madenci (2010)](https://doi.org/10.1016/j.tafmec.2010.08.001).


## Export settings

The results of a simulation are exported as `vtu` and `jld2` files, to make post-processing with ParaView and Julia easy.
With [`ExportSettings`](@ref), the file path and the frequency of the export can be set.


## Predefined cracks

Predefined Cracks are defined with [`PreCrack`](@ref).


## Boundary and initial conditions

Currently, 4 different conditions can be set:

*Boundary conditions*
1. [`VelocityBC`](@ref): velocity boundary conditions
2. [`ForceDensityBC`](@ref): force density boundary conditions
3. [`PosDepVelBC`](@ref): position dependend velocity boundary conditions

*Initial conditions*
1. [`VelocityIC`](@ref): velocity boundary conditions