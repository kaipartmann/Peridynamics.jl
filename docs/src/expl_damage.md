# [Damage in peridynamics formulations](@id expl_dmg)

For damage simulations with peridynamics, some slight extensions to the material formulations need to be made.

## Bond-based formulation

In bond-based peridynamics, the pairwise force function is expanded by the factor $d^{ij}$, which states whether the bond between points $i$ and $j$ is intact ($d^{ij}=1$) or damaged ($d^{ij}=0$).
Therefore the pairwise force function reads
```math
 \boldsymbol{f} = d^{ij} \, c \, \varepsilon^{ij} \, \boldsymbol{n} 
```
with the bond failure quantity $d^{ij} \in \{0;1\}$.

## State-based formulation

In all state-based peridynamic formulations, damage is introduced by using a failure considering influence function:
```math
 \omega_d = d^{ij} \, \omega
```
with the bond failure quantity $d^{ij} \in \{0;1\}$.


## Continuum-kinematics-inspired formulation

For continuum-kinematics-inspired peridynamics there are three different bond failure factors, one for each kind of interaction.
For one-neighbor interactions it is similar to the failure in the bond-based formulation. Here the internal force density due to one-neighbor interactions is
```math
    \boldsymbol{b}_1^{{\mathrm{int}},i} = \int_{\mathcal{H}_1^i} d^{ij} C_1 \left( \frac{l^{ij}}{L^{ij}} - 1 \right) \frac{\boldsymbol{\Delta x}^{ij}}{l^{ij}} \; \mathrm{d} V_1^i
```
with the bond failure quantity $d^{ij} \in \{0;1\}$.

In two-neighbor interactions a factor describing the failure of the considered area element is included:
```math
    \boldsymbol{b}_{2}^{\mathrm{int}, \, i} = 
2 \, C_2 \int_{\mathcal{H}_2^i} d^{ijk} \left( \frac{a^{ijk}}{A^{ijk}} - 1 \right)
\frac{\boldsymbol{\Delta x}^{ik} \times \boldsymbol{a}^{ijk}}{a^{ijk}} \; \mathrm{d} V_2^i 
```
with the two-neighbor interaction failure quantity $d^{ijk} \in \{0;1\}$.

For three-neighbor interactions the internal force density eventually reads
```math
\boldsymbol{b}_{3}^{\mathrm{int}, \, i} = 
3 \, C_3 \int_{\mathcal{H}_3^i} d^{ijkl} \left( \frac{\left|{v^{ijkl}}\right|}{\left|{V^{ijkl}}\right|} - 1 \right)
\frac{\left(\boldsymbol{\Delta x}^{ik} \times \boldsymbol{\Delta x}^{il}\right) v^{ijkl}}{\left|{v^{ijkl}}\right|} \; \mathrm{d} V_3^i
```
with the three-neighbor interaction failure quantity $d^{ijkl} \in \{0;1\}$.

| Size | Symbol |      Unit |
|:--------|:-------------|:------------|
| Bond failure quantity |      $d^{ij} \in \{0;1\}$      | $[-]$ |
| Two-neighbor interaction failure quantity |      $d^{ijk} \in \{0;1\}$      | $[-]$ |
| Three-neighbor interaction failure quantity |      $d^{ijkl} \in \{0;1\}$      | $[-]$ |
