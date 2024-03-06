# [Continuum-kinematics-inspired](@id continuumbased)

Continuum-kinematics-inspired peridynamics (CPD) is a formulation that is supposed to deliver more freedom in specifying material parameters.
The internal force density is calculated as the sum of three types of point interactions which are one-, two- and three-neighbor interactions [Javili2019](@cite):

```math
\boldsymbol{b}^{\mathrm{int},i} = \boldsymbol{b}_1^{\mathrm{int},i} + \boldsymbol{b}_2^{\mathrm{int},i} + \boldsymbol{b}_3^{\mathrm{int},i}
``` 

## One-neighbor interactions

```@raw html
<img src="/assets/PD_OneNI_1.png?raw=true" width="250"/>
```

One-neighbor interactions in CPD correspond to the bonds in [bond-based](@ref bondbased) peridynamics, but there is a slightly different way to calculate the internal forces.

First, the neighborhood volume is determined:
```math
V_\mathcal{H}^i = \beta^i \cdot \frac {4} {3} \cdot\pi\cdot\delta^3
```

Here $\beta^i\in [0,1]$ is a factor for the completeness of the neighborhood that takes incomplete point families at the surface into account (see figure 1).

```@raw html
<img src="/assets/OberflaechenEffekt_3.png?raw=true" width="250"/>
```

Now, the effective one-neighbor volume can be calculated
```math
V_1^i = \frac{V_\mathcal{H}^i}{N_1^i}
```
with the number of interactions $N_1^i$.

The internal force density is determined by
```math
    \boldsymbol{b}_1^{{\mathrm{int}},i} = \int_{\mathcal{H}_1^i} d^{ij} C_1 \left( \frac{l^{ij}}{L^{ij}} - 1 \right) \frac{\boldsymbol{\Delta x}^{ij}}{l^{ij}} \; \mathrm{d} V_1^i
```
with the parameters:

| size | symbol |      unit |
|:--------|:-------------|:------------|
| bond failure |      $d^{ij} \in \{0;1\}$      | $[-]$ |
| material constant | $C_1$ |      $[\frac{\mathrm{kg}}{\mathrm{m}^5\mathrm{s}^2}]$ |
| relative length measures | $L^{ij}$, $l^{ij}$ | $[\mathrm{m}]$ |
| effective one-neighbor volume   |  $V_1^i = \frac{V_\mathcal{H}^i}{N_1^i}$                 | $[\mathrm{m}^3]$ |
| neighborhood volume | $V_\mathcal{H}^i = \beta^i \cdot \frac {4}{3} \cdot\pi\cdot\delta^3$ | $[\mathrm{m}^3]$ |
| neighborhood completeness   |     $\beta^i\in [0,1]$      | $[-]$ |
| number of one-neighbor interactions   |      $N_1^i$      | $[-]$ |


## Two-neighbor interactions

```@raw html
<img src="/assets/PD_TwoNI_1.png?raw=true" width="250"/>
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
2 \, C_2 \int_{\mathcal{H}_2^i} d^{ijk} \left( \frac{a^{ijk}}{A^{ijk}} - 1 \right)
\frac{\boldsymbol{\Delta x}^{ik} \times \boldsymbol{a}^{ijk}}{a^{ijk}} \; \mathrm{d} V_2^i \; .
```

| size | symbol |      unit |
|:--------|:-------------|:------------|
| material constant | $C_2$ |      $[\frac{\mathrm{kg}}{\mathrm{m}^9\mathrm{s}^2}]$ |
| (bond) failure |      $d^{ijk} \in \{0;1\}$      | $[-]$ |
| relative area measures | $A^{ijk}$, $a^{ijk}$, $\boldsymbol{a}^{ijk}$ | $[\mathrm{m}^2]$ |
| effective two-neighbor volume   |  $V_2^i = \frac{\left(V_\mathcal{H}^i\right)^2}{N_2^i}$  | $[\mathrm{m}^6]$ |
| number of two-neighbor interactions   |      $N_2^i$      | $[-]$ |

## Three-neighbor interactions

```@raw html
<img src="/assets/PD_ThreeNI_1.png?raw=true" width="250"/>
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
3 \, C_3 \int_{\mathcal{H}_3^i} d^{ijkl} \left( \frac{\left|{v^{ijkl}}\right|}{\left|{V^{ijkl}}\right|} - 1 \right)
\frac{\left(\boldsymbol{\Delta x}^{ik} \times \boldsymbol{\Delta x}^{il}\right) v^{ijkl}}{\left|{v^{ijkl}}\right|} \; \mathrm{d} V_3^i
```

with the material constant $C_3$ is used.

| size | symbol |      unit |
|:--------|:-------------|:------------|
| material constant | $C_3$ |      $[\frac{\mathrm{kg}}{\mathrm{m}^{13}\mathrm{s}^2}]$ |
| (bond) failure |      $d^{ijkl} \in \{0;1\}$      | $[-]$ |
| relative volume measures | $V^{ijkl}$, $v^{ijkl}$ | $[\mathrm{m}^3]$ |
| effective three-neighbor volume |  $V_3^i = \frac{ \left(V_\mathcal{H}^i\right)^3}{N_3^i}$ | $[\mathrm{m}^9]$ |
| number of three-neighbor interactions |      $N_3^i$      | $[-]$ |
