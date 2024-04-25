# [Ordinary State-based](@id statebased)

The state-based peridynamics formulation considers not only the deformation of the bonds of one material point as the [bond-based model](@ref bondbased) does, but also the states of all neighbors to calculate the internal force density $\boldsymbol{b}^{\mathrm{int},i}$ [Silling2007](@cite):

```math
\boldsymbol{b}^{\mathrm{int}} (\boldsymbol{X}^i,t) = \int_{\mathcal{H}_i} \boldsymbol{t}^i - \boldsymbol{t}^j \; \mathrm{d}V^j \; .
```
with the force vector states $\boldsymbol{t}^i=\boldsymbol{t}(\boldsymbol{\Delta X}^{ij}, t)$ and $\boldsymbol{t}^j=\boldsymbol{t}(-\boldsymbol{\Delta X}^{ij}, t)$, which characterize the state of each bond at time $t$.

To determine the force vector state, the extension scalar state $e$ is required first.
It describes the change in length of a bond due to deformation, consisting of a spherical and a deviatoric component:
```math
e \langle \boldsymbol{\Delta X}^{ij} \rangle = |\boldsymbol{\Delta x}^{ij}|-|\boldsymbol{\Delta X}^{ij}|
= e^{sph} \langle \boldsymbol{\Delta X}^{ij} \rangle + e^{dev} \langle \boldsymbol{\Delta X}^{ij} \rangle \; .
```

Next, the weighted volume $m_i$ is calculated as
```math
m_i = \left(\omega \langle \boldsymbol{\Delta X}^{ij} \rangle \cdot | \boldsymbol{\Delta X}^{ij} |\right) \bullet | \boldsymbol{\Delta X}^{ij} | = \int_{\mathcal{H}_i} \omega \langle \boldsymbol{\Delta X}^{ij} \rangle \cdot | \boldsymbol{\Delta X}^{ij} |^2 \mathrm{d} V_j \; .
```
Here $\omega$ is an influence function that gives greater weight to shorter bonds. [Silling2007](@cite)

Then the dilatation $\theta_i$ is needed, which is determined with the weighted volume $m_i$:
```math
\theta_i = \frac{3}{m_i} \left(\omega \langle \boldsymbol{\Delta X}^{ij} \rangle \cdot | \boldsymbol{\Delta X}^{ij} |\right) \bullet e \langle \boldsymbol{\Delta X}^{ij} \rangle
= \frac{3}{m_i} \int_{\mathcal{H}_i} \omega \langle \boldsymbol{\Delta X}^{ij} \rangle \cdot | \boldsymbol{\Delta X}^{ij} | \cdot e \langle \boldsymbol{\Delta X}^{ij} \rangle \mathrm{d} V_j \; .
```

The spherical component of the extension scalar state $e$ results in 
```math
e^{sph} \langle \boldsymbol{\Delta X}^{ij} \rangle = \frac{\theta_i \cdot |\boldsymbol{\Delta X}^{ij}|}{3}
```
and the deviatoric component in
```math
e^{dev} \langle \boldsymbol{\Delta X}^{ij} \rangle = |\boldsymbol{\Delta x}^{ij}|-|\boldsymbol{\Delta X}^{ij}| - \frac{\theta_i \cdot |\boldsymbol{\Delta X}^{ij}|}{3} \; .
```

With the previously determined variables, the force vector state $\boldsymbol{t}$ can now be calculated: [Silling2007](@cite)
```math
\underline{t}_i \boldsymbol{t}^i \langle \boldsymbol{\Delta X}^{ij} \rangle = \frac{\kappa \cdot \theta_i}{m_i} \cdot \omega \langle \boldsymbol{\Delta X}^{ij} \rangle \cdot | \boldsymbol{\Delta X}^{ij} | + \frac{15 \cdot \mu}{m_i} \cdot \omega \langle \boldsymbol{\Delta X}^{ij} \rangle \cdot e^{dev} \langle \boldsymbol{\Delta X}^{ij} \rangle \; .
```
