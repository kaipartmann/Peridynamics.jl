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