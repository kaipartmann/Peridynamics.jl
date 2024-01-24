# Non-Ordinary State-based

Non-ordinary state based formulations have been developed to extend state based peridynamics.
Hereafter, the correspondence formulation of non-ordinary state based peridynamics is considered, which uses an elastic model from the classical theory. [Silling2007](@cite)

First, the symmetric shape tensor is calculated:
```math
\boldsymbol{K} = \boldsymbol{K}(\boldsymbol{X}) = \int_\mathcal{H} \omega \, \boldsymbol{\Delta X} \otimes \boldsymbol{\Delta X} \; \mathrm{d}V' \; .
```
Here, $\omega$ is an influence function to weigh points differently.
The deformation gradient is thus approximated as [Silling2007](@cite)
```math
\boldsymbol{F} = \boldsymbol{F}(\boldsymbol{X},t) = \left(\int_\mathcal{H} \omega \, \boldsymbol{\Delta x} \otimes \boldsymbol{\Delta X} \; \mathrm{d}V'\right) \boldsymbol{K}^{-1} \; .
```

Using the deformation gradient, now the first Piola Kirchhoff stress tensor can be determined with the Helmholtz energy density $\Psi$:
```math
\boldsymbol{P} = \boldsymbol{P}(\boldsymbol{X},t) = \frac{\partial \Psi}{\partial \boldsymbol{F}} \; .%= \boldsymbol{F} \boldsymbol{S} \; .
```

Using the calculated variables, the force vector state can now be determined by [Silling2007](@cite)
```math
\boldsymbol{t} = \omega \boldsymbol{P}  \boldsymbol{K}^{-1} \boldsymbol{\Delta X} \; .
```

