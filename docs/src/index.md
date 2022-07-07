```@raw html
<img src="https://github.com/kfrb/Peridynamics.jl/blob/main/docs/src/assets/logo.png?raw=true" width="450" />
```
# Peridynamics

Documentation for the [Peridynamics](https://github.com/kfrb/Peridynamics.jl) julia package.

## Overview

Peridynamics is a non-local formulation of continuum mechanics.
Material points map the continuum, and the relative displacements and forces are described by an integral equation that is also fulfilled for discontinuities. Therefore, peridynamics is ideally suited for dynamic fracture simulations with many cracks!

This package aims at making it easy for everyone to perform peridynamics simulations with a high-level API and the speed of Julia.
Currently, only a limited feature set is supported, but many more are in the pipeline, so stay tuned for upcoming versions!

#### Features
- Bond-based peridynamics, see [Silling (2000)](https://doi.org/10.1016/S0022-5096(99)00029-0)
- Import and convert meshes from Abaqus (Supported element types: `Tet4`, `Hex8`)
- Explicit time integration with Velocity Verlet algorithm
- Adaptive dynamic relaxation for quasistatic analysis, see [Kilic and Madenci (2010)](https://doi.org/10.1016/j.tafmec.2010.08.001)
- Multi-body contact analysis with short-range forces, see [Silling and Askari (2005)](https://doi.org/10.1016/j.compstruc.2004.11.026)

#### Incomplete list of upcoming features
- Volume and surface correction
- Ordinary and non-ordinary state-based peridynamics, see [Silling et al. (2007)](https://link.springer.com/article/10.1007/s10659-007-9125-1)
- Continuum-kinematics-based peridynamics, see [Javili, McBride and Steinmann (2019)](https://doi.org/10.1016/j.jmps.2019.06.016)

## Installation

To install, use Julia's built-in package manager. Open the Julia REPL and type `]` to enter the package mode and install Peridynamics as follows:

```julia-repl
pkg> add Peridynamics
```

## Getting Started

We recommend looking at the [manual](@ref manual) and the [tutorials section](@ref TensileTest.jl) to start working with this package!

Have fun! 😃