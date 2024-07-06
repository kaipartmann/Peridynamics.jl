<img src="docs/src/assets/logo.png" width="450" />

# Peridynamics.jl
A high-level Julia package for parallel peridynamics simulations

**Documentation:**\
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kaipartmann.github.io/Peridynamics.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kaipartmann.github.io/Peridynamics.jl/dev/)

**Build status:**\
[![Build Status](https://github.com/kaipartmann/Peridynamics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kaipartmann/Peridynamics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kaipartmann/Peridynamics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kaipartmann/Peridynamics.jl)

**Code quality:**\
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2)](https://github.com/SciML/SciMLStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

**Citation:**\
[![DOI](https://zenodo.org/badge/503281781.svg)](https://zenodo.org/badge/latestdoi/503281781)

## Main Features
- 🎯 Dynamic and quasi-static analysis
- 🪨 Multiple peridynamics formulations and material models
- 🎳 Multibody contact simulations
- 🤓 User friendly API that captures many errors before submitting simulations
- 🚀 Enhanced HPC capabilities with MPI and multiple threads

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

## Tutorials

<table align="center" border="0">
  <tr>
  </tr>
  <tr>
    <td align="center">
      <a href="https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_tension_static/">
        <figcaption>Tensile test quasi-static</figcaption><img src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/ac69d8aa-526d-436a-aa0c-820a1f42bcca" style="width: 90% !important;"/><br>
      </a>
    </td>
    <td align="center">
      <a href="https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_tension_dynfrac/">
        <figcaption>Tensile test dynamic</figcaption><img src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/dda2b7b3-d44b-41a9-b133-6d1b548df1c1" style="width: 90% !important;"/><br>
      </a>
    </td>
  </tr>
  <tr>
  </tr>
  <tr>
    <td align="center">
      <a href="https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_tension_precrack/">
        <figcaption>Tension with predefined crack</figcaption><img src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/9f627d2d-44b5-43a3-94cd-9d34894fd142" style="width: 90% !important;"/><br>
      </a>
    </td>
    <td align="center">
      <a href="https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_logo/">
        <figcaption>Peridynamics.jl logo</figcaption><img src="docs/src/assets/logo.png" style="width: 90% !important;"/><br>
      </a>
    </td>
  </tr>
</table>

## Authors

- <a href="https://orcid.org/0000-0002-5238-4355">Kai Partmann (University of Siegen) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
- <a href="https://orcid.org/0009-0004-9195-0112">Manuel Dienst (University of Siegen) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
- <a href="https://orcid.org/0000-0002-2213-8401">Kerstin Weinberg (University of Siegen) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>

## Acknowledgements
<img src=https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/0d14a65b-4e05-4408-8107-59ac9c1477d2 width=500>

The authors gratefully acknowledge the support of the Deutsche Forschungsgemeinschaft (DFG) under the project WE2525-14/1.

The authors gratefully acknowledge the computing time provided to them on the high-performance computer Noctua 2 at the NHR Center PC2. These are funded by the Federal Ministry of Education and Research and the state governments participating on the basis of the resolutions of the GWK for the national highperformance computing at universities (www.nhr-verein.de/unsere-partner).
