# Peridynamics

A high-level Julia package for multithreaded peridynamics simulations

## What is peridynamics?
Peridynamics is a non-local formulation of continuum mechanics where material points represent the continuum, and the relative displacements and forces are governed by an integro-differential equation that allows discontinuities. As a result, peridynamics is particularly well-suited for dynamic fracture simulations involving numerous cracks.

## Who can benefit from this package?
This package is designed for anyone interested in performing peridynamics simulations. While the current feature set provides a solid foundation, we are continuously working to enhance and expand the capabilities of `Peridynamics.jl`. We encourage you to open an issue or submit a pull request to share your feedback or contribute to making this package even more valuable to the community!

## Installation

To install `Peridynamics.jl`, follow these steps:

1. Install Julia from the [official Julia website](https://julialang.org/) if you haven't already.

2. Launch Julia and open the Julia REPL.

3. Enter the package manager by pressing `]` in the REPL.

4. In the package manager, type:
   ```
   add Peridynamics
   ```

5. Press `Backspace` or `Ctrl + C` to exit the package manager.

## How-to guides
#### How to perform a
- [Single body analysis](@ref howto_single_body_analysis)
- [Contact analysis](@ref howto_contact_analysis)

#### How to define
- [Point clouds](@ref howto_pointclouds)
- [Predefined cracks](@ref howto_precracks)
- [Material formulations](@ref howto_matformulations)

#### How to
- [Visualize results with ParaView](@ref visualization)

## Tutorials

```@raw html
<div class="tutorial-grid"> 
```

```@raw html
<div class="tutorial-element"> 
```

### [Mode I tension quasi-static](@ref tutorial_tension_static)
[![](assets/tension_static.gif)](@ref tutorial_tension_static)

```@raw html
</div> 
```

```@raw html
<div class="tutorial-element"> 
```

### [Mode I tension dynamic](@ref tutorial_tension_dynfrac)
[![](assets/tension_dynfrac.gif)](@ref tutorial_tension_dynfrac)

```@raw html
</div>
```

```@raw html
<div class="tutorial-element"> 
```

### [Mode I tension dynamic with predefined crack](@ref tutorial_tension_precrack)
[![](assets/tension_precrack_damage.gif)](@ref tutorial_tension_precrack)

```@raw html
</div> 
```

```@raw html
<div class="tutorial-element"> 
```

### [Peridynamics.jl logo](@ref tutorial_logo)
[![](assets/logo.gif)](@ref tutorial_logo)

```@raw html
</div>
```
```@raw html
</div>
```

## Authors

```@raw html
<ul>
<li><a href="https://orcid.org/0000-0002-5238-4355">Kai Partmann (University of Siegen) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a></li>
<li><a href="https://orcid.org/0009-0004-9195-0112">Manuel Dienst (University of Siegen) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a></li>
<li><a href="https://orcid.org/0000-0002-2213-8401">Kerstin Weinberg (University of Siegen) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a></li>
</ul>
```

## Acknowledgement
```@raw html
<img src=https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/0d14a65b-4e05-4408-8107-59ac9c1477d2 width=500>
```
The authors gratefully acknowledge the support of the Deutsche Forschungsgemeinschaft (DFG) under the project WE2525-14/1.

The authors gratefully acknowledge the computing time provided to them on the high-performance computer Noctua 2 at the NHR Center PC2. These are funded by the Federal Ministry of Education and Research and the state governments participating on the basis of the resolutions of the GWK for the national highperformance computing at universities (www.nhr-verein.de/unsere-partner).
