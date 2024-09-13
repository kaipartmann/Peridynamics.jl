# [Ordinary state-based peridynamics](@id expl_osb)

The state-based peridynamics formulation considers not only the deformation of the bonds of one material point, but also the states of all neighbors to calculate the internal force density $\boldsymbol{b}^{\mathrm{int},i}$ as
```math
\boldsymbol{b}^{\mathrm{int},i} = \boldsymbol{b}^{\mathrm{int}} (\boldsymbol{X}^i,t) = \int_{\mathcal{H}_i} \boldsymbol{t}^i - \boldsymbol{t}^j \; \mathrm{d}V^j \; ,
```
with the force vector states $\boldsymbol{t}^i=\boldsymbol{t}(\boldsymbol{\Delta X}^{ij}, t)$ and $\boldsymbol{t}^j=\boldsymbol{t}(-\boldsymbol{\Delta X}^{ij}, t)$, which characterize the state of each bond at time $t$ [Silling2007](@cite).
To determine the force vector states, the weighted volume $m_i$ is calculated first as
```math
m_i = m \left( \boldsymbol{X}^i \right) = \int_{\mathcal{H}_i} \omega \, | \boldsymbol{\Delta X}^{ij} |^2 \, \mathrm{d} V^j \; .
```
Here, $\omega$ is the influence function that gives a greater influence to neighbors near the root point. [Silling2007](@cite)
Then the dilatation $\theta_i$ is needed, which is defined with the weighted volume $m_i$ as
```math
\theta_i = \theta \left( \boldsymbol{X}^i \right) = \frac{3}{m_i} \int_{\mathcal{H}_i} \omega \, | \boldsymbol{\Delta X}^{ij} | \, \left( |\boldsymbol{\Delta x}^{ij}|-|\boldsymbol{\Delta X}^{ij}| \right) \mathrm{d} V^j \; .
```
With the previously determined variables, the force vector state $\boldsymbol{t}^i$ is defined as 
```math
\boldsymbol{t}^i \left( \boldsymbol{\Delta X}^{ij} \right) = \frac{K \, \theta_i}{m_i} \, \omega \, | \boldsymbol{\Delta X}^{ij} | + \frac{15 \, G}{m_i} \, \omega \, \left( |\boldsymbol{\Delta x}^{ij}|-|\boldsymbol{\Delta X}^{ij}| - \frac{\theta_i \, |\boldsymbol{\Delta X}^{ij}|}{3} \right) \; .
```
with shear modulus $G$ and bulk modulus $K$ [Silling2007](@cite).

| Size | Symbol |      Unit |
|:--------|:-------------|:------------|
| Internal force density | $\boldsymbol{b}^{\mathrm{int},i}$ | $\left[\frac{\mathrm{kg}}{\mathrm{m}^2\mathrm{s}^2}\right]$ |
| Force vector state | $ \boldsymbol{t}^i $ | $\left[\frac{\mathrm{kg}}{\mathrm{m}^5\mathrm{s}^2}\right]$ |
| Volume of point $j$ | $V^j$ | $\left[\mathrm{m}^3\right]$
| Bond in $\mathcal{B}_0$ |      $\boldsymbol{\Delta X}^{ij}$     | $[\mathrm{m}]$ |
| Bond in $\mathcal{B}_t$ |      $\boldsymbol{\Delta x}^{ij}$     | $[\mathrm{m}]$ |
| Bond length in $\mathcal{B}_0$ |      $\left\|\boldsymbol{\Delta X}^{ij}\right\|$     | $[\mathrm{m}]$ |
| Bond length in $\mathcal{B}_t$ |      $\left\|\boldsymbol{\Delta x}^{ij}\right\|$     | $[\mathrm{m}]$ |
| Influence function | $\omega$ |      $[-]$ |
| Weighted volume |  $ m_i $ | $\left[\mathrm{m}^5\right]$ |
| Dilatation | $\theta_i$ |      $[-]$ |
| Shear modulus |      $G$      | $\left[\frac{\mathrm{kg}}{\mathrm{m}\mathrm{s}^2}\right]$ |
| Bulk modulus |      $K$      | $\left[\frac{\mathrm{kg}}{\mathrm{m}\mathrm{s}^2}\right]$ |
