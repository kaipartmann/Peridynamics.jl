# Public API

Here, all exported types and methods are listed, which the user needs when using the package.

```@meta

```

```@index
Pages = ["api_public.md"]
```

## Types
```@docs
BBMaterial
OSBMaterial
NOSBMaterial
CKIMaterial
Body
MultibodySetup
VelocityVerlet
Job
```

## Functions
```@docs
point_set!
point_sets
n_points
failure_permit!
material!
velocity_bc!
velocity_ic!
forcedensity_bc!
precrack!
contact!
read_vtk
uniform_box
submit
```