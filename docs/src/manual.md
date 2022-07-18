# [Manual](@id manual)

!!! tip "Typical workflow"
    `Peridynamics.jl` simulations follow this short workflow:

    1. Create a job
    2. Submit the job
   

## How to create a job

Currently two different simulation job types are available:

- [`PDSingleBodyAnalysis`](@ref)
  - Define the spatial discretization with [`PointCloud`](@ref)
  - Define the material, e.g. [`BondBasedMaterial`](@ref)
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

### Bond-based
Currently only the [`BondBasedMaterial`](@ref) is implemented.
More material models will be added in the future!

#### Example:
- Horizon $\delta = 0.008\,\mathrm{m}$
- Density $\rho = 2700\,\mathrm{kg}\,\mathrm{m}^{-3}$
- Youngs modulus $E = 70 \times 10^9 \, \mathrm{Pa}$
- Griffith's parameter $G_c = 100 \, \mathrm{N} \, \mathrm{m}^{-1}$
```julia
material = BondBasedMaterial(; horizon=0.008, rho=2700, E=70e9, Gc=100)
```


## Temporal discretization

The temporal discretization is controlled with the [`TimeDiscretization`](@ref) type.
By default, explicit time integration with velocity verlet algorithm is used.
Quasistatic simulations can be performed using the adaptive dynamic relaxation algorithm by [Kilic and Madenci (2010)](https://doi.org/10.1016/j.tafmec.2010.08.001).

#### Example 1:
Velocity verlet algorithm with 500 time steps, the step size is calculated.
```julia
TimeDiscretization(500)
```

#### Example 2:
Velocity verlet algorithm with 500 time steps and a step size of $\Delta t = 5 \times 10^{-7} \, \mathrm{s}$.
```julia
TimeDiscretization(500, 5e-7)
```

## Export settings

The simulation results are exported as `vtu` and `jld2` files to make post-processing with ParaView and Julia straightforward.
With [`ExportSettings`](@ref), the file path and the export frequency can be set.

#### Example 1:
Export into the directory `"results"` every 10'th time step.
```julia
ExportSettings("results", 10)
```

#### Example 2:
Export nothing.
```julia
ExportSettings()
```

## Predefined cracks

Predefined cracks are defined with [`PreCrack`](@ref).
For an example, see [`CrackedPlateUnderTension.jl`](@ref cracked-plate-under-tension).


## Boundary & initial conditions

For more information, see the documentation for each type:

- [`VelocityBC`](@ref): velocity boundary conditions
- [`ForceDensityBC`](@ref): force density boundary conditions
- [`PosDepVelBC`](@ref): position dependend velocity boundary conditions
- [`VelocityIC`](@ref): velocity initial conditions