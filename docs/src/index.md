# Peridynamics

A high-level Julia package for multithreaded peridynamics simulations

| **Cite the package** |**Downloads counter**|
|:---|:---|
| ... | [![Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/Peridynamics/label:-color:blue)](https://pkgs.genieframework.com?packages=Peridynamics) |

## What is peridynamics?
Peridynamics is a non-local formulation of continuum mechanics where material points represent the continuum, and the relative displacements and forces are governed by an integro-differential equation that allows discontinuities. As a result, peridynamics is particularly well-suited for dynamic fracture simulations involving numerous cracks.

## Who can benefit from this package?
This package is designed for anyone interested in performing peridynamics simulations. With its intuitive high-level API and the incredible performance of Julia, users can seamlessly delve into the world of peridynamics. While the current feature set provides a solid foundation, we are continuously working to enhance and expand the capabilities of `Peridynamics.jl`. We encourage you to open an issue or submit a pull request to share your feedback and make this package even more valuable to the community!

## How-to guides

- [Spatial discretization](@ref spatial_discretization)
- [Material models](@ref mat_models)
- [Temporal discretization](@ref time_discretization)
- [Conditions](@ref conditions)
- [Jobs / PD-Analysis](@ref jobs)
- [Workflow](@ref workflow)
- [Visualization with ParaView](@ref visualization)
- [Read VTK results](@ref vtk_reader)

## Tutorials

```@raw html
<div class="tutorial-grid"> 
```

```@raw html
<div class="tutorial-element"> 
```

### [Tensile test](@ref tensile_test)
[![](assets/TensileTest.png)](@ref tensile_test)

```@raw html
</div> 
```

```@raw html
<div class="tutorial-element"> 
```

### [Cracked plate under tension](@ref cracked-plate-under-tension)
[![](assets/CrackedPlateUnderTension2000.png)](@ref cracked-plate-under-tension)

```@raw html
</div> 
```

```@raw html
<div class="tutorial-element"> 
```

### [Peridynamics.jl logo](@ref logo)
[![](assets/logo.gif)](@ref logo)

```@raw html
</div> 
```

```@raw html
</div> 
```

## Authors

- Kai Partmann (University of Siegen)
- Kerstin Weinberg (University of Siegen)

## Acknowledgement
```@raw html
<img src=https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/0d14a65b-4e05-4408-8107-59ac9c1477d2 width=500>
```
The authors gratefully acknowledge the support of the Deutsche Forschungsgemeinschaft (DFG) under the project WE2525-14/1.
