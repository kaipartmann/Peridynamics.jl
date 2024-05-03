# [Ordinary State-based](@id statebased)

The state-based peridynamics formulation considers not only the deformation of the bonds of one material point as the [bond-based model](@ref bondbased) does, but also the states of all neighbors to calculate the internal force density $\boldsymbol{b}^{\mathrm{int},i}$ [Silling2007](@cite):

```math
\boldsymbol{b}^{\mathrm{int},i} = \boldsymbol{b}^{\mathrm{int}} (\boldsymbol{X}^i,t) = \int_{\mathcal{H}_i} \boldsymbol{t}^i - \boldsymbol{t}^j \; \mathrm{d}V^j \; .
```
with the force vector states $\boldsymbol{t}^i=\boldsymbol{t}(\boldsymbol{\Delta X}^{ij}, t)$ and $\boldsymbol{t}^j=\boldsymbol{t}(-\boldsymbol{\Delta X}^{ij}, t)$, which characterize the state of each bond at time $t$.

To determine the force vector states, the weighted volume $m_i$ is calculated first as
```math
m_i = m \left( \boldsymbol{X}^i \right) = \int_{\mathcal{H}_i} \omega \cdot | \boldsymbol{\Delta X}^{ij} |^2 \, \mathrm{d} V^j \; .
```
Here $\omega$ is an influence function that gives greater weight to shorter bonds. [Silling2007](@cite)

Then the dilatation $\theta_i$ is needed, which is determined with the weighted volume $m_i$:
```math
\theta_i = \theta \left( \boldsymbol{X}^i \right) = \frac{3}{m_i} \int_{\mathcal{H}_i} \omega \cdot | \boldsymbol{\Delta X}^{ij} | \cdot \left( |\boldsymbol{\Delta x}^{ij}|-|\boldsymbol{\Delta X}^{ij}| \right) \mathrm{d} V^j \; .
```

With the previously determined variables, the force vector state $\boldsymbol{t}^i$ can now be calculated: [Silling2007](@cite)
```math
\boldsymbol{t}^i \left( \boldsymbol{\Delta X}^{ij} \right) = \frac{\kappa \cdot \theta_i}{m_i} \cdot \omega \cdot | \boldsymbol{\Delta X}^{ij} | + \frac{15 \cdot \mu}{m_i} \cdot \omega \cdot \left( |\boldsymbol{\Delta x}^{ij}|-|\boldsymbol{\Delta X}^{ij}| - \frac{\theta_i \cdot |\boldsymbol{\Delta X}^{ij}|}{3} \right) \; .
```
