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

```julia
submit(job)
```
