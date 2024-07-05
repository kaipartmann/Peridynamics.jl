# [Bond-based](@id bondbased)

The initial version of peridynamics is the bond-based (BB) formulation. [Silling2000](@cite)

Here, a pairwise force function $\boldsymbol{f}$ is defined and calculated for each bond of two material points, which depends on the strain of the bond and is aligned in its direction: 

```math
 \boldsymbol{f} = d^{ij} \cdot c \cdot \varepsilon^{ij} \cdot \boldsymbol{n} \; .
```

with

| size | symbol |      unit |
|:--------|:-------------|:------------|
| bond failure |      $d^{ij} \in \{0;1\}$      | $[-]$ |
| micro-modulus constant [Silling2005](@cite) |  $c = \frac{18 \cdot \kappa}{\pi \cdot \delta^4}$ | $[\frac{\mathrm{kg}}{\mathrm{m}^5\mathrm{s}^2}]$ |
| bond strain | $\varepsilon^{ij} = \frac{l^{ij}-L^{ij}}{L^{ij}}$ |      $[-]$ |
| bond length in $\mathcal{B}_0$ |      $L^{ij} =\left\|\boldsymbol{\Delta X}^{ij}\right\|$     | $[\mathrm{m}]$ |
| bond length in $\mathcal{B}_t$ |      $l^{ij} =\left\|\boldsymbol{\Delta x}^{ij}\right\|$     | $[\mathrm{m}]$ |
| direction vector |      $\boldsymbol{n} = \frac{\boldsymbol{\Delta x}^{ij}}{l^{ij}}$      | $[-]$ |

In this expression, bond failure $d^{ij}$ represents whether the bond between points $i$ and $j$ is intact ($d^{ij}=1$)
or damaged ($d^{ij}=0$). The direction vector $\boldsymbol{n}$ is oriented in the direction of the bond. 
Furthermore, the micro-modulus constant [Silling2005](@cite)
```math
c = \frac{18 \cdot \kappa}{\pi \cdot \delta^4} \;
```
and the strain of the bond [Silling2005a](@cite)
```math
\varepsilon^{ij} = \frac{l^{ij}-L^{ij}}{L^{ij}}
``` 
with bond lengths $L^{ij} =\left|\boldsymbol{\Delta X}^{ij}\right|$ and $l^{ij} =\left|\boldsymbol{\Delta x}^{ij}\right|$ are used.

To get the resulting body forces, now the force function is integrated over the whole body:

```math
\boldsymbol{b}^{\mathrm{int},i} = \boldsymbol{b}^{\mathrm{int}} (\boldsymbol{X} ^ {i} , t) = \int_{\mathcal{H}_i} \boldsymbol{f} \; \mathrm{d}V^j \; .
```

