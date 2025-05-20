# Public API

```@meta

```

```@contents
Pages = ["public_api_reference.md"]
```

## Material models
```@docs
BBMaterial
DHBBMaterial
OSBMaterial
CMaterial
CRMaterial
BACMaterial
CKIMaterial
RKCMaterial
RKCRMaterial
```

## System or material related types
```@docs
CriticalStretch
NoCorrection
EnergySurfaceCorrection
ZEMSilling
ZEMWan
LinearElastic
NeoHooke
MooneyRivlin
SaintVenantKirchhoff
const_one_kernel
linear_kernel
cubic_b_spline_kernel
```

## Discretization
```@docs
Body
MultibodySetup
point_set!
point_sets
no_failure!
material!
velocity_bc!
velocity_ic!
forcedensity_bc!
precrack!
contact!
uniform_box
uniform_sphere
uniform_cylinder
round_sphere
round_cylinder
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