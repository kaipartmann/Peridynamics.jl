# [Continuum-kinematics-inspired peridynamics](@id expl_cki)

Continuum-kinematics-inspired peridynamics (CPD) is a formulation that is supposed to deliver more freedom in specifying material parameters.
The internal force density is calculated as the sum of three types of point interactions which are one-, two- and three-neighbor interactions [Javili2019](@cite):

```math
\boldsymbol{b}^{\mathrm{int},i} = \boldsymbol{b}_1^{\mathrm{int},i} + \boldsymbol{b}_2^{\mathrm{int},i} + \boldsymbol{b}_3^{\mathrm{int},i} \; .
``` 

| Size | Symbol |      Unit |
|:--------|:-------------|:------------|
| Internal force density | $\boldsymbol{b}^{\mathrm{int},i}$ | $\left[\frac{\mathrm{kg}}{\mathrm{m}^2\mathrm{s}^2}\right]$ |
| Force density shares due to one-, two- & three-neighbor interactions | $\boldsymbol{b}_1^{\mathrm{int},i}$ , $\boldsymbol{b}_2^{\mathrm{int},i}$ , $\boldsymbol{b}_3^{\mathrm{int},i}$ | $\left[\frac{\mathrm{kg}}{\mathrm{m}^2\mathrm{s}^2}\right]$ |

## One-neighbor interactions

```@raw html
<img src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/0380b1b8-4527-4f38-a20e-435a6f1c8ba1" width="250"/>
```

One-neighbor interactions in CPD correspond to the bonds in [bond-based](@ref expl_bb) peridynamics, but there is a slightly different way to calculate the internal forces.

First, the neighborhood volume is determined:
```math
V_\mathcal{H}^i = \beta^i \, \frac {4} {3} \, \pi \, \delta^3 \; .
```

Here $\beta^i\in [0,1]$ is a factor for the completeness of the neighborhood that takes incomplete point families at the surface into account (see figure 1).

```@raw html
<img src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/b899c8d3-e358-4d4d-b52f-13b6c0af747b" width="350"/>
```

Now, the effective one-neighbor volume can be calculated
```math
V_1^i = \frac{V_\mathcal{H}^i}{N_1^i}
```
with the number of interactions $N_1^i$.

The internal force density is determined by
```math
    \boldsymbol{b}_1^{{\mathrm{int}},i} = \int_{\mathcal{H}_1^i} C_1 \left( \frac{l^{ij}}{L^{ij}} - 1 \right) \frac{\boldsymbol{\Delta x}^{ij}}{l^{ij}} \; \mathrm{d} V_1^i
```
with the parameters:

| Size | Symbol |      Unit |
|:--------|:-------------|:------------|
| Neighborhood volume | $V_\mathcal{H}^i$ | $[\mathrm{m}^3]$ |
| Neighborhood completeness   |     $\beta^i\in [0,1]$      | $[-]$ |
| Effective one-neighbor volume   |  $V_1^i$                 | $[\mathrm{m}^3]$ |
| Number of one-neighbor interactions   |      $N_1^i$      | $[-]$ |
| Material constant | $C_1$ |      $[\frac{\mathrm{kg}}{\mathrm{m}^5\mathrm{s}^2}]$ |
| Relative length measures | $L^{ij}$, $l^{ij}$ | $[\mathrm{m}]$ |


## Two-neighbor interactions

```@raw html
<img src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/5b340634-7c3f-4ddf-a056-06c31564077c" width="250"/>
```

For two-neighbor interactions, the deformation of the area spanned by point $i$ and two of its neighbors $j$ and $k$ is analyzed to calculate the internal force density. For this, relative area measures are defined:

```math
    A^{ijk}=\left| \boldsymbol{\Delta X}^{ij} \times \boldsymbol{\Delta X}^{ik} \right| \; , \qquad a^{ijk}=\left| \boldsymbol{\Delta x}^{ij} \times \boldsymbol{\Delta x}^{ik} \right| \; , \qquad \boldsymbol{a}^{ijk}= \boldsymbol{\Delta x}^{ij} \times \boldsymbol{\Delta x}^{ik} \; .
```

Other sizes needed to identify the force density are the material constant $C_2$ and the effective two-neighbor volume
```math
    V_2^i = \frac{\left(V_\mathcal{H}^i\right)^2}{N_2^i}
```
with the number of interactions $N_2$. 

The internal force density induced by two-neighbor interactions is 

```math
    \boldsymbol{b}_{2}^{\mathrm{int}, \, i} = 
2 \, C_2 \int_{\mathcal{H}_2^i} \left( \frac{a^{ijk}}{A^{ijk}} - 1 \right)
\frac{\boldsymbol{\Delta x}^{ik} \times \boldsymbol{a}^{ijk}}{a^{ijk}} \; \mathrm{d} V_2^i \; .
```

| Size | Symbol |      Unit |
|:--------|:-------------|:------------|
| Relative area measures | $A^{ijk}$, $a^{ijk}$, $\boldsymbol{a}^{ijk}$ | $[\mathrm{m}^2]$ |
| Effective two-neighbor volume   |  $V_2^i$  | $[\mathrm{m}^6]$ |
| Number of two-neighbor interactions   |      $N_2^i$      | $[-]$ |
| Material constant | $C_2$ |      $[\frac{\mathrm{kg}}{\mathrm{m}^9\mathrm{s}^2}]$ |

## Three-neighbor interactions

```@raw html
<img src="https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/e908c804-a6d4-4bf6-b75a-988f33989213" width="250"/>
```

Three-neighbor interactions regard the volume defined by the bond vectors between point $i$ and its three neighbors $j$, $k$ and $l$:

```math
V^{ijkl} = \left(\boldsymbol{\Delta X}^{ij} \times \boldsymbol{\Delta X}^{ik}\right) \cdot \boldsymbol{\Delta X}^{il}  \;,\qquad
    v^{ijkl} = \left(\boldsymbol{\Delta x}^{ij} \times \boldsymbol{\Delta x}^{ik}\right) \cdot \boldsymbol{\Delta x}^{il}  \;.
```
Additionally, the effective three-neighbor volume
```math
    V_3^i = \frac{ \left(V_\mathcal{H}^i\right)^3}{N_3^i} \; .
```
is defined.
For the internal force density of three-neighbor interactions, the equation

```math
\boldsymbol{b}_{3}^{\mathrm{int}, \, i} = 
3 \, C_3 \int_{\mathcal{H}_3^i} \left( \frac{\left|{v^{ijkl}}\right|}{\left|{V^{ijkl}}\right|} - 1 \right)
\frac{\left(\boldsymbol{\Delta x}^{ik} \times \boldsymbol{\Delta x}^{il}\right) v^{ijkl}}{\left|{v^{ijkl}}\right|} \; \mathrm{d} V_3^i
```

with the material constant $C_3$ is used.

| Size | Symbol |      Unit |
|:--------|:-------------|:------------|
| Relative volume measures | $V^{ijkl}$, $v^{ijkl}$ | $[\mathrm{m}^3]$ |
| Effective three-neighbor volume |  $V_3^i$ | $[\mathrm{m}^9]$ |
| Number of three-neighbor interactions |      $N_3^i$      | $[-]$ |
| Material constant | $C_3$ |      $[\frac{\mathrm{kg}}{\mathrm{m}^{13}\mathrm{s}^2}]$ |
