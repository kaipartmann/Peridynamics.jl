```@raw html
<p align="center">
  <img src="assets/logo.png" width="300" />
  <br>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/817c7bd4-9c02-4cc4-ac66-998c0f5e95e2">
    <source media="(prefers-color-scheme: light)" srcset="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/70c24007-5aa9-460f-9a97-c67b1df32750">
    <img alt="The Peridynamics.jl logo" src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/70c24007-5aa9-460f-9a97-c67b1df32750" width="400">
  </picture>
</p>
```

A high-level Julia package for parallel peridynamics simulations

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
- [Simulations with MPI](@ref)
- [Visualize results with ParaView](@ref visualization)

## Tutorials

```@raw html
<div class="tutorial-grid">
```

```@raw html
<div class="tutorial-element">
   <a href="https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_tension_static/">
      <figcaption>Tensile test quasi-static</figcaption><br><img src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/ac69d8aa-526d-436a-aa0c-820a1f42bcca" style="width: 90% !important;"/>
   </a>
</div>
```

```@raw html
<div class="tutorial-element">
   <a href="https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_tension_dynfrac/">
      <figcaption>Tensile test dynamic</figcaption><br><img src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/dda2b7b3-d44b-41a9-b133-6d1b548df1c1" style="width: 90% !important;"/>
   </a>
</div>
```


```@raw html
<div class="tutorial-element">
   <a href="https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_tension_precrack/">
      <figcaption>Tension with predefined crack</figcaption><br><img src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/9f627d2d-44b5-43a3-94cd-9d34894fd142" style="width: 90% !important;"/>
   </a>
</div>
```

```@raw html
<div class="tutorial-element">
   <a href="https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_logo/">
      <figcaption>The old logo</figcaption><br><img src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/5439e112-9088-49a3-bb01-aff541adc0f8" style="width: 90% !important;"/>
   </a>
</div>
```

```@raw html
<div class="tutorial-element">
   <a href="https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_kalthoff-winkler_dynfrac/">
      <figcaption>Kalthoff Winkler</figcaption><br><img src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/6dc362ef-4997-4327-9bc1-41350fac2dc1" style="width: 90% !important;"/>
   </a>
</div>
```

```@raw html
<div class="tutorial-element">
   <a href="https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_cylinder/">
      <figcaption>Fragmenting cylinder</figcaption><br><img src="https://github.com/user-attachments/assets/58e11123-6143-4e13-8642-7e30c9e6c86d" style="width: 90% !important;"/>
   </a>
</div>
```

```@raw html
<div class="tutorial-element">
   <a href="https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_wave_in_bar/">
      <figcaption>Wave propagation</figcaption><br><img src="https://github.com/user-attachments/assets/7fa65fd4-38d8-46cb-833f-990417211d17" style="width: 90% !important;"/>
   </a>
</div>
```

```@raw html
<div class="tutorial-element">
   <a href="https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_wave_interface/">
      <figcaption>Wave propagation across interface</figcaption><br><img src="https://github.com/user-attachments/assets/082f635f-caf2-40db-938e-e4a98e2f3915" style="width: 90% !important;"/>
   </a>
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

## Acknowledgements
```@raw html
<img src=https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/0d14a65b-4e05-4408-8107-59ac9c1477d2 width=500>
```
The authors gratefully acknowledge the support of the Deutsche Forschungsgemeinschaft (DFG) under the project WE2525-14/1.

The support of Carsten Bauer and Xin Wu from PC2 with the design of the internal structure regarding parallel performance is gratefully acknowledged.

The authors gratefully acknowledge the computing time provided to them on the high-performance computer Noctua 2 at the NHR Center PC2. These are funded by the Federal Ministry of Education and Research and the state governments participating on the basis of the resolutions of the GWK for the national highperformance computing at universities (www.nhr-verein.de/unsere-partner).
