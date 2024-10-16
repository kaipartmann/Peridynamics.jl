<p align="center">
  <img src="docs/src/assets/logo.png" width="300" />
  <br>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/817c7bd4-9c02-4cc4-ac66-998c0f5e95e2">
    <source media="(prefers-color-scheme: light)" srcset="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/70c24007-5aa9-460f-9a97-c67b1df32750">
    <img alt="The Peridynamics.jl logo" src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/70c24007-5aa9-460f-9a97-c67b1df32750" width="400">
  </picture>
</p>

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
[![DOI](https://proceedings.juliacon.org/papers/10.21105/jcon.00165/status.svg)](https://doi.org/10.21105/jcon.00165)
[![DOI](https://zenodo.org/badge/503281781.svg)](https://zenodo.org/badge/latestdoi/503281781)

## Main Features
- 🎯 Dynamic and quasi-static analysis
- 🪨 Multiple peridynamics formulations and material models
- 🎳 Multibody contact simulations
- 🤓 User friendly API that captures many errors before submitting simulations
- 🚀 Enhanced HPC capabilities with MPI or multithreading

## Installation

`Peridynamics.jl` is a registered Julia package, so you can install it by just typing
```
add Peridynamics
```
in the julia package manager. Please take a look at the [documentation](https://kaipartmann.github.io/Peridynamics.jl/stable/index#Installation) for more details on the installation.

## Tutorials

<p align="center">
  <img src="docs/src/assets/tutorial_overview.svg" width="100%" height="900" />
</p>

## Usage
To run the dynamic tensile test simulation shown above, just 7 lines of code are needed:
```julia
body = Body(BBMaterial(), "TensileTestMesh.inp")
material!(body; horizon=0.01, rho=2700, E=70e9, Gc=100)
velocity_bc!(t -> -0.6, body, :bottom, 1)
velocity_bc!(t -> 0.6, body, :top, 1)
vv = VelocityVerlet(steps=500)
job = Job(body, vv; path="results/tension_dynamic")
submit(job)
```
Take a look at the [tutorial of the tensile test](https://kaipartmann.github.io/Peridynamics.jl/stable/generated/tutorial_tension_dynfrac/) for more details on this example.

If you want to run this example with multithreading, just start Julia with more than 1 thread.
To use MPI, you can create a script containing the **same code without changes** and run it with:
```bash
mpiexec -n 6 julia --project path/to/script.jl
```
Please take a look at the [how-to guide on MPI](https://kaipartmann.github.io/Peridynamics.jl/dev/howto_mpi/) for more details.

## Authors

- <a href="https://orcid.org/0000-0002-5238-4355">Kai Partmann (University of Siegen) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
- <a href="https://orcid.org/0009-0004-9195-0112">Manuel Dienst (University of Siegen) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
- <a href="https://orcid.org/0000-0002-2213-8401">Kerstin Weinberg (University of Siegen) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>

## Acknowledgements
<img src=https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/0d14a65b-4e05-4408-8107-59ac9c1477d2 width=500>

The authors gratefully acknowledge the support of the Deutsche Forschungsgemeinschaft (DFG) under the project WE2525-14/1.

The support of Carsten Bauer and Xin Wu from PC2 with the design of the internal structure regarding parallel performance is gratefully acknowledged.

The authors gratefully acknowledge the computing time provided to them on the high-performance computer Noctua 2 at the NHR Center PC2. These are funded by the Federal Ministry of Education and Research and the state governments participating on the basis of the resolutions of the GWK for the national highperformance computing at universities (www.nhr-verein.de/unsere-partner).
