# [Bond-based peridynamics](@id expl_bb)

The initial version of peridynamics is the bond-based (BB) formulation. [Silling2000](@cite)

Here, a pairwise force function $\boldsymbol{f}$ is defined and calculated for each bond of two material points, which depends on the strain of the bond and is aligned in its direction: 

```math
 \boldsymbol{f} = c \, \varepsilon^{ij} \, \boldsymbol{n} \; .
```

Here the micro-modulus constant [Silling2005](@cite)
```math
c = \frac{18 \, \kappa}{\pi \, \delta^4} \;
```
and the strain of the bond [Silling2005a](@cite)
```math
\varepsilon^{ij} = \frac{l^{ij}-L^{ij}}{L^{ij}}
``` 
with bond lengths $L^{ij} =\left|\boldsymbol{\Delta X}^{ij}\right|$ and $l^{ij} =\left|\boldsymbol{\Delta x}^{ij}\right|$ are used.

The direction vector
```math
\boldsymbol{n} = \frac{\boldsymbol{\Delta x}^{ij}}{l^{ij}}
``` 
is oriented in the direction of the bond. 

To get the resulting body forces, now the force function is integrated over the whole body:

```math
\boldsymbol{b}^{\mathrm{int},i} = \boldsymbol{b}^{\mathrm{int}} (\boldsymbol{X} ^ {i} , t) = \int_{\mathcal{H}_i} \boldsymbol{f} \; \mathrm{d}V^j \; .
```

| Size | Symbol |      Unit |
|:--------|:-------------|:------------|
| Pairwise force function | $\boldsymbol{f}$ | $\left[\frac{\mathrm{kg}}{\mathrm{m}^5\mathrm{s}^2}\right]$ |
| Micro-modulus constant [Silling2005](@cite) |  $c$ | $\left[\frac{\mathrm{kg}}{\mathrm{m}^5\mathrm{s}^2}\right]$ |
| Bond strain | $\varepsilon^{ij}$ |      $[-]$ |
| Bond in $\mathcal{B}_0$ |      $\boldsymbol{\Delta X}^{ij}$     | $[\mathrm{m}]$ |
| Bond in $\mathcal{B}_t$ |      $\boldsymbol{\Delta x}^{ij}$     | $[\mathrm{m}]$ |
| Bond length in $\mathcal{B}_0$ |      $L^{ij}$     | $[\mathrm{m}]$ |
| Bond length in $\mathcal{B}_t$ |      $l^{ij}$     | $[\mathrm{m}]$ |
| Direction vector |      $\boldsymbol{n}$      | $[-]$ |
| Volume of point $j$ | $V^j$ | $\left[\mathrm{m}^3\right]$
| Internal force density | $\boldsymbol{b}^{\mathrm{int},i}$ | $\left[\frac{\mathrm{kg}}{\mathrm{m}^2\mathrm{s}^2}\right]$ |
