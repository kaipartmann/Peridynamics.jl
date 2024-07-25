# [Damage in peridynamics formulations](@id expl_dmg)

## Bond-based formulation

```math
 \boldsymbol{f} = d^{ij} \, c \, \varepsilon^{ij} \, \boldsymbol{n} \; .
```

with

| size | symbol |      unit |
|:--------|:-------------|:------------|
| bond failure |      $d^{ij} \in \{0;1\}$      | $[-]$ |

In this expression, bond failure $d^{ij}$ represents whether the bond between points $i$ and $j$ is intact ($d^{ij}=1$)
or damaged ($d^{ij}=0$).

## Ordinary state-based formulation

## Non-ordinary state-based formulation

## Continuum-kinematics-inspired formulation

```math
    \boldsymbol{b}_1^{{\mathrm{int}},i} = \int_{\mathcal{H}_1^i} d^{ij} C_1 \left( \frac{l^{ij}}{L^{ij}} - 1 \right) \frac{\boldsymbol{\Delta x}^{ij}}{l^{ij}} \; \mathrm{d} V_1^i
```
with the parameters:

| size | symbol |      unit |
|:--------|:-------------|:------------|
| bond failure |      $d^{ij} \in \{0;1\}$      | $[-]$ |

```math
    \boldsymbol{b}_{2}^{\mathrm{int}, \, i} = 
2 \, C_2 \int_{\mathcal{H}_2^i} d^{ijk} \left( \frac{a^{ijk}}{A^{ijk}} - 1 \right)
\frac{\boldsymbol{\Delta x}^{ik} \times \boldsymbol{a}^{ijk}}{a^{ijk}} \; \mathrm{d} V_2^i \; .
```

| size | symbol |      unit |
|:--------|:-------------|:------------|
| material constant | $C_2$ |      $[\frac{\mathrm{kg}}{\mathrm{m}^9\mathrm{s}^2}]$ |
| (area element) failure |      $d^{ijk} \in \{0;1\}$      | $[-]$ |

```math
\boldsymbol{b}_{3}^{\mathrm{int}, \, i} = 
3 \, C_3 \int_{\mathcal{H}_3^i} d^{ijkl} \left( \frac{\left|{v^{ijkl}}\right|}{\left|{V^{ijkl}}\right|} - 1 \right)
\frac{\left(\boldsymbol{\Delta x}^{ik} \times \boldsymbol{\Delta x}^{il}\right) v^{ijkl}}{\left|{v^{ijkl}}\right|} \; \mathrm{d} V_3^i
```

with the material constant $C_3$ is used.

| size | symbol |      unit |
|:--------|:-------------|:------------|
| material constant | $C_3$ |      $[\frac{\mathrm{kg}}{\mathrm{m}^{13}\mathrm{s}^2}]$ |
| (volume element) failure |      $d^{ijkl} \in \{0;1\}$      | $[-]$ |

```math
\boldsymbol{b}^{\mathrm{int},i} = \boldsymbol{b}_1^{\mathrm{int},i} + \boldsymbol{b}_2^{\mathrm{int},i} + \boldsymbol{b}_3^{\mathrm{int},i} \; .
``` 
