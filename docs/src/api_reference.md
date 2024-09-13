# Peridynamics.jl public API

```@meta

```

```@contents
Pages = ["api_reference.md"]
```

## Material models
```@docs
BBMaterial
OSBMaterial
NOSBMaterial
CKIMaterial
```

## System related types
```@docs
NoCorrection
EnergySurfaceCorrection
```

## Discretization
```@docs
Body
MultibodySetup
point_set!
point_sets
failure_permit!
material!
velocity_bc!
velocity_ic!
forcedensity_bc!
precrack!
contact!
uniform_box
uniform_sphere
n_points
```

## Preprocessing & simulation setup
```@docs
read_inp
mpi_isroot
force_mpi_run!
force_threads_run!
enable_mpi_timers!
disable_mpi_timers!
enable_mpi_progress_bars!
reset_mpi_progress_bars!
@mpitime
@mpiroot
```

## Solving
```@docs
VelocityVerlet
DynamicRelaxation
Job
submit
```

## Postprocessing
```@docs
read_vtk
process_each_export
```