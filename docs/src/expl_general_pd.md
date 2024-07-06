# [Basics of peridynamics theory](@id expl_pd_basics)

Peridynamics is a nonlocal continuum mechanics formulation, which was introduced by Silling [Silling2000](@cite).
It has gained increased popularity as an approach for modeling fracture.
The deformation of the solid is described by integro-differential equations that are also fulfilled for discontinuities, making it very capable of modeling crack propagation and fragmentation with large displacements.
Much peridynamics research has been done in recent years, summarized in various review papers and books [Diehl2019,Javili2019Review,Madenci2014](@cite).

```@raw html
<img src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/728da1f0-4750-4ab6-a430-9b206e475577" width="300"/>
```

Typically, in peridynamics the continuum is discretized by material points.
Points interact only with other points inside of their specified *neighborhood* or *point family* $\mathcal{H}$, which is defined as the set of points inside a sphere with the radius $\delta$, also named the *horizon*.
The interaction of the point $\boldsymbol{X}$ with its neighbor $\boldsymbol{X}'$ is called *bond* and defined as
```math
\boldsymbol{\Delta X} = \boldsymbol{X}' - \boldsymbol{X} \; .
```
The equation of motion reads
```math
\varrho \, \boldsymbol{\ddot{u}}(\boldsymbol{X},t) = \boldsymbol{b}^{\mathrm{int}}(\boldsymbol{X},t) + \boldsymbol{b}^{\mathrm{ext}}(\boldsymbol{X},t) \; ,
```
with the mass density $\varrho$, the point acceleration vector $\boldsymbol{\ddot{u}}$, and the point force density vectors $\boldsymbol{b}^{\mathrm{int}}$ and $\boldsymbol{b}^{\mathrm{ext}}$.
Various material formulations of peridynamics exist for the calculation of the internal force density $\boldsymbol{b}^{\mathrm{int}}$, and all of them are based on the nonlocal interactions between material points.

The general internal force density for state-based peridynamics is defined as
```math
\boldsymbol{b}^{\mathrm{int}} (\boldsymbol{X},t) = \int_\mathcal{H} \boldsymbol{t} - \boldsymbol{t}' \; \mathrm{d}V' \; ,
```
with the *force vector states* $\boldsymbol{t}=\boldsymbol{t}(\boldsymbol{\Delta X}, t)$ and $\boldsymbol{t}'=\boldsymbol{t}(-\boldsymbol{\Delta X}, t)$.
In the first original [bond-based formulation of peridynamics](@ref expl_bb), the force vector states $\boldsymbol{t}$ and  $\boldsymbol{t}'$ have the same value and opposite direction.
This implies intrinsic limitation to only one material parameter and in consequence to restrictions on the Poisson's ratio [Silling2007,Trageser2020](@cite).
To overcome these restrictions, state-based peridynamics was established.
In the [ordinary state-based peridynamics](@ref expl_osb), the deformation states of neighboring points also influence the internal force density [Silling2007](@cite).
This leads to force vector states which are still collinear but not of same value anymore.

Further developments are summarized as non-ordinary state-based peridynamics.
A recent development in this regard is [continuum-kinematics-inspired peridynamics](@ref expl_cki) [Javili2019](@cite).
Another peridynamic formulation is the [local continuum consistent correspondence formulation](@ref expl_nosb) of non-ordinary state-based peridynamics, where an elastic model from the classical local material theory can be used to calculate the internal force density.
