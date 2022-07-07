# TensileTest.jl

## Quickstart

##### Code:
```julia
using Peridynamics
pointcloud = read_inp("TensileTestMesh.inp")
material = BondBasedMaterial(; horizon=0.008, rho=2700, E=70e9, Gc=100)
boundary_conditions = [
    VelocityBC(t -> -0.8, pointcloud.point_sets["bottom"], 1),
    VelocityBC(t -> 0.8, pointcloud.point_sets["top"], 1),
]
job = PDSingleBodyAnalysis(;
    name="TensileTest",
    pc=pointcloud,
    mat=material,
    bcs=boundary_conditions,
    td=TimeDiscretization(500),
    es=ExportSettings("results/TensileTest", 10),
)
submit(job)
```

##### Results:
```@raw html
<img src="https://github.com/kfrb/Peridynamics.jl/blob/main/docs/src/assets/TensileTest.png?raw=true" width="600" />
```
(Visualization made with [ParaView](https://www.paraview.org))

## Step-by-step
The script for this tutorial can be found here: [TensileTest.jl](https://github.com/kfrb/Peridynamics.jl/blob/main/examples/TensileTest.jl).

First, we have to load the `Peridynamics.jl` package.

```julia
using Peridynamics
```
Now we create a `PointCloud` by converting the mesh from Abaqus. [(Download the mesh)](https://github.com/kfrb/Peridynamics.jl/blob/main/examples/models/TensileTestMesh.inp)
```julia
pointcloud = read_inp("TensileTestMesh.inp")
```
```@raw html
<img src="https://github.com/kfrb/Peridynamics.jl/blob/main/docs/src/assets/TensileTestMesh.png?raw=true" width="600" />
```

Next, we need to specify material parameters. We define:
- Horizon $\delta = 0.008\,\mathrm{m}$
- Density $\rho = 2700\,\mathrm{kg}\,\mathrm{m}^{-3}$
- Youngs modulus $E = 70 \times 10^9 \, \mathrm{Pa}$
- Griffith's parameter $G_c = 100 \, \mathrm{N} \, \mathrm{m}^{-1}$
```julia
material = BondBasedMaterial(; horizon=0.008, rho=2700, E=70e9, Gc=100)
```
As loading condition for the specimen, a constant velocity of $0.8 \, \mathrm{m}\,\mathrm{s}^{-1}$ in $x$-direction is set for the bottom and top.
Note, that element sets defined in Abaqus are converted to `point_sets` of the `PointCloud`.
```julia
boundary_conditions = [
    VelocityBC(t -> -0.8, pointcloud.point_sets["bottom"], 1),
    VelocityBC(t -> 0.8, pointcloud.point_sets["top"], 1),
]
```
With these at hand, a `PDSingleBodyAnalysis` can be constructed and submitted for calculation.
We set the number of time steps for the explicit time integration to 500 with `TimeDiscretization(500)` and define, that the results of our calculation should be saved in the directory `"results/TensileTest"` every 10'th time step.
```julia
job = PDSingleBodyAnalysis(;
    name="TensileTest",
    pc=pointcloud,
    mat=material,
    bcs=boundary_conditions,
    td=TimeDiscretization(500),
    es=ExportSettings("results/TensileTest", 10),
)
submit(job)
```

The output of `submit(job)` looks as follows:
```
======================================================================
PERIDYNAMIC SIMULATION ON 8 THREADS
======================================================================
Peridynamic single body analysis: TensileTest
Material parameters:
  - Number of material points [-]:                               16900
  - Material model:                                  BondBasedMaterial
  - Horizon δ [m]:                                               0.008
  - Density ρ [kg/m³]:                                            2700
  - Young's modulus E [N/m²]:                                    7e+10
  - Poisson ratio ν [-]:                                          0.25
  - Critical bond stretch εc[-]:                           0.000385758
Interactions:
  - Number of bonds [-]:                                       6667240
Total memory used by body [MB]:                                273.856
Export setup:
  - Export frequency:                                               10
  - Export file name:         examples/results/TensileTest/TensileTest
Time discretization:
  - Time step Δt [s]:                                      5.13584e-07
  - Number of time steps [-]:                                      500
  - Simulation time horizon [s]:                           0.000256792
Time integration... 100%|██████████████████████████████| Time: 0:00:13
✓ Simulation completed after 16.7321 seconds
Results:
  - Max. abs. x-displacement [m]:                          0.000434543
  - Max. abs. y-displacement [m]:                          0.000326501
  - Max. abs. z-displacement [m]:                          0.000230936
  - Max. damage [-]:                                                 1
```

## No-fail-zone

As an improvement to the previous results, we can additionally set a no-fail-zone.
It is a common problem for bond-based peridynamics, that to much damage occurs next to the points where boundary conditions apply.
To disable failure of the points in the bottom and the top of the specimen, add these two lines after the definition of the `pointcloud`:
```julia
pointcloud.failure_flag[pointcloud.point_sets["bottom"]] .= false
pointcloud.failure_flag[pointcloud.point_sets["top"]] .= false
```
As a result, we have reduced the damage in the bottom and top of the specimen:
```@raw html
<img src="https://github.com/kfrb/Peridynamics.jl/blob/main/docs/src/assets/TensileTestNoFailZone.png?raw=true" width="600" />
```