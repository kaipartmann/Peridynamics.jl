<img src="docs/src/assets/logo.png" width="450" />

# Peridynamics.jl
A high-level Julia package for multithreaded peridynamics simulations

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

## Getting Started

We recommend looking at the [how-to guides](https://kaipartmann.github.io/Peridynamics.jl/stable/) and the [tutorials](https://kaipartmann.github.io/Peridynamics.jl/stable/) to start working with this package!

## Tutorials Overview

<style>
.tutorial-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    grid-gap: 20px;
    padding: 25px;
}

.tutorial-element {
    text-align: center;
    /* background-color: #f9f9f9; */
    background-color: rgb(240,240,240);;
    padding: 20px;
    border-radius: 10px;
    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
    transition: transform 0.3s, box-shadow 0.3s;
}

.tutorial-element:hover {
    transform: scale(1.05);
    box-shadow: 0 8px 16px rgba(0, 0, 0, 0.2);
}

.tutorial-element img {
    width: 100%;
    max-width: 300px;
    height: auto;
    border-radius: 10px;
    margin-bottom: 10px;
}

.tutorial-element a {
    color: #333;
    text-decoration: none;
    font-weight: bold;
}

.tutorial-element a:hover {
    color: #007bff;
}
</style>

<div class="tutorial-grid">

<div class="tutorial-element">
<a href="link_to_your_documentation#tutorial_tension_static">Tensile test quasi-static</a><br>
<a href="link_to_your_documentation#tutorial_tension_static">
    <img src="docs/src/assets/tension_static.gif" alt="Tensile test quasi-static">
</a>
</div>

<div class="tutorial-element">
<a href="link_to_your_documentation#tutorial_tension_dynfrac">Tensile test dynamic</a><br>
<a href="link_to_your_documentation#tutorial_tension_dynfrac">
    <img src="docs/src/assets/tension_dynfrac.gif" alt="Tensile test dynamic">
</a>
</div>

<div class="tutorial-element">
<a href="link_to_your_documentation#tutorial_tension_precrack">Tension with predefined crack</a><br>
<a href="link_to_your_documentation#tutorial_tension_precrack">
    <img src="docs/src/assets/tension_precrack_damage.gif" alt="Tension with predefined crack">
</a>
</div>

<div class="tutorial-element">
<a href="link_to_your_documentation#tutorial_logo">Peridynamics.jl logo</a><br>
<a href="link_to_your_documentation#tutorial_logo">
    <img src="docs/src/assets/logo.gif" alt="Peridynamics.jl logo">
</a>
</div>

</div>


## Authors

- <a href="https://orcid.org/0000-0002-5238-4355">Kai Partmann (University of Siegen) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
- <a href="https://orcid.org/0009-0004-9195-0112">Manuel Dienst (University of Siegen) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
- <a href="https://orcid.org/0000-0002-2213-8401">Kerstin Weinberg (University of Siegen) <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>

## Acknowledgement
<img src=https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/0d14a65b-4e05-4408-8107-59ac9c1477d2 width=500>

The authors gratefully acknowledge the support of the Deutsche Forschungsgemeinschaft (DFG) under the project WE2525-14/1.

The authors gratefully acknowledge the computing time provided to them on the high-performance computer Noctua 2 at the NHR Center PC2. These are funded by the Federal Ministry of Education and Research and the state governments participating on the basis of the resolutions of the GWK for the national highperformance computing at universities (www.nhr-verein.de/unsere-partner).
