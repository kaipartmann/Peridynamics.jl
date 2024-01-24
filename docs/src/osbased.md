# Ordinary State-based

The state-based peridynamics formulation considers not only the deformation of the bonds of one material point, but also the states of all neighbors, to calculate the internal force density $\boldsymbol{b}^{\mathrm{int},i}$ [Silling2007](@cite):

```math
\boldsymbol{b}^{\mathrm{int}} (\boldsymbol{X},t) = \int_\mathcal{H} \boldsymbol{t} - \boldsymbol{t}' \; \mathrm{d}V' 
```
with the force vector states $\boldsymbol{t}=\boldsymbol{t}(\boldsymbol{\Delta X}, t)$ and $\boldsymbol{t}'=\boldsymbol{t}(-\boldsymbol{\Delta X}, t)$, which characterize the state of each bond at time $t$.

To determine the force vector state, the extension scalar state $\underline{e}$ is required first.
It describes the change in length of a bond due to deformation, consisting of a spherical and a deviatoric component:
```math
\underline{e} \langle \boldsymbol{\Delta X} \rangle \left( = |\underline{\boldsymbol{Y}}_k \langle \boldsymbol{\Delta X} \rangle | - | \boldsymbol{\underline{X}}_k \langle \boldsymbol{\Delta X} \rangle | \right) = |\boldsymbol{\Delta x}|-|\boldsymbol{\Delta X}|
= \underline{e}^{sph} \langle \boldsymbol{\Delta X} \rangle + \underline{e}^{dev} \langle \boldsymbol{\Delta X} \rangle
```

Next, the weighted volume $m_k$ is calculated as
```math
m_k = (\underline{w} \langle \boldsymbol{\Delta X} \rangle \cdot | \boldsymbol{\Delta X} |) \bullet | \boldsymbol{\Delta X} | = \int_{\mathcal{H}_k} \underline{w} \langle \boldsymbol{\Delta X} \rangle \cdot | \boldsymbol{\Delta X} |^2 \mathrm{d} V_j \; .
```
Here $\omega$ is an influence function that gives greater weight to shorter bonds. [Silling2007](@cite)

Then the dilataion $\theta_k$ is needed, which is determined with the weighted volume $m_k$:
```math
\theta_k = \frac{3}{m_k} (\underline{w} \langle \boldsymbol{\Delta X} \rangle \cdot | \boldsymbol{\Delta X} |) \bullet \underline{e} \langle \boldsymbol{\Delta X} \rangle
= \frac{3}{m_k} \int_{\mathcal{H}_k} \underline{w} \langle \boldsymbol{\Delta X} \rangle \cdot | \boldsymbol{\Delta X} | \cdot \underline{e} \langle \boldsymbol{\Delta X} \rangle \mathrm{d} V_j
```

The spherical component of the extension scalar state $\underline{e}$ results in 
```math
\underline{e}^{sph} \langle \boldsymbol{\Delta X} \rangle = \frac{\theta_k \cdot |\boldsymbol{\Delta X}|}{3}
```
and the deviatoric component in
```math
\underline{e}^{dev} \langle \boldsymbol{\Delta X} \rangle = |\boldsymbol{\Delta x}|-|\boldsymbol{\Delta X}| - \frac{\theta_k \cdot |\boldsymbol{\Delta X}|}{3} \; .
```

With the previously determined variables, the force vector state $\boldsymbol{t}$ can now be calculated:
```math
\underline{t}_k \langle \boldsymbol{\Delta X} \rangle = \frac{\kappa \cdot \theta_k}{m_k} \cdot \underline{w} \langle \boldsymbol{\Delta X} \rangle \cdot | \boldsymbol{\Delta X} | + \frac{15 \cdot \mu}{m_k} \cdot \underline{w} \langle \boldsymbol{\Delta X} \rangle \cdot \underline{e}^{dev} \langle \boldsymbol{\Delta X} \rangle
```
