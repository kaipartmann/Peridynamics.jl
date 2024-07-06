# [Non-ordinary state-based peridynamics](@id expl_nosb)

Non-ordinary state based-formulations have been developed to extend [state-based peridynamics](@ref expl_osb).
Hereafter, the correspondence formulation of non-ordinary state based peridynamics is considered, which uses an elastic model from the classical theory. [Silling2007](@cite)

First, the symmetric shape tensor is calculated:
```math
\boldsymbol{K}^i = \boldsymbol{K}(\boldsymbol{X}^i) = \int_{\mathcal{H}_i} \omega \, \boldsymbol{\Delta X}^{ij} \otimes \boldsymbol{\Delta X}^{ij} \; \mathrm{d}V^j \; .
```
Here, $\omega$ is an influence function to weigh points differently.
The deformation gradient is thus approximated as [Silling2007](@cite)
```math
\boldsymbol{F}^i = \boldsymbol{F}(\boldsymbol{X}^i,t) = \left(\int_{\mathcal{H}_i} \omega \, \boldsymbol{\Delta x}^{ij} \otimes \boldsymbol{\Delta X}^{ij} \; \mathrm{d}V^j\right) \left(\boldsymbol{K}^i\right)^{-1} \; .
```

Using the deformation gradient, now the first Piola Kirchhoff stress tensor can be determined with the Helmholtz energy density $\Psi$:
```math
\boldsymbol{P}^i = \boldsymbol{P}(\boldsymbol{X}^i,t) = \frac{\partial \Psi}{\partial \boldsymbol{F}^i} \; .%= \boldsymbol{F} \boldsymbol{S} \; .
```

Using the calculated variables, the force vector state can now be determined by [Silling2007](@cite)
```math
\boldsymbol{t}^i = \omega \boldsymbol{P}^i  \left(\boldsymbol{K}^i\right)^{-1} \boldsymbol{\Delta X}^{ij} \; .
```
