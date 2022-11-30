# Public API

Here, all exported types and methods are listed, which the user needs when using the package.

```@meta
CurrentModule = Peridynamics
```

```@index
Pages = ["public_api.md"]
```

## Types
```@docs
PointCloud
BondBasedMaterial
PreCrack
TimeDiscretization
ExportSettings
BodySetup
Contact
PDSingleBodyAnalysis
PDContactAnalysis
VelocityBC
VelocityIC
PosDepVelBC
ForceDensityBC
SimResult
```

## Functions
```@docs
calc_stable_user_timestep
read_inp
read_vtk
sphere_radius
submit
```