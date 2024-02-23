using Peridynamics #hide

# # [Material formulations](@id howto_matformulations)
# In peridynamics multiple material formulations exist. In this package, currently two
# peridynamics formulations are implemented:
#   - [`BBMaterial`](@ref)
#   - [`CPDMaterial`](@ref)
# Please refer to
#   - [Silling (2000)](https://doi.org/10.1016/S0022-5096(99)00029-0)
#   - [Silling et al. (2007)](https://link.springer.com/article/10.1007/s10659-007-9125-1)
#   - [Javili, McBride and Steinmann (2019)](https://doi.org/10.1016/j.jmps.2019.06.016)
# for further details on the theoretical aspects of these peridynamics formulations.

# ## Bond-based peridynamics
# The bond-based formulation was the first peridynamics model. It has the built-in
# limitation to a Poisson-ration of $\nu=1/4$. Otherwise you can specify the horizon
#  $\delta$, density $\rho$, Young's modulus $E$ and either the critical energy release rate
# $\mathcal{G}_c$ or the critical bond strain $\varepsilon_c$. All other needed parameters
# are then calculated.

# ## Examples bond-based material
#-
mat = BBMaterial(horizon=1, rho=8e-6, E=2.1e5, Gc=2)

#-
mat = BBMaterial(horizon=1, rho=8e-6, E=2.1e5, epsilon_c=0.01)


# ## Continuum-kinematics-based peridynamics
# The continuum-kinematics-based (also called continuum-kinematics-inspired) formulation
# has in it's core formulas three material parameters $C_1$, $C_2$, and $C_3$.
# [Ekiz et al.](https://doi.org/10.1016/j.mechmat.2022.104417) found, that
# ```math
# C_1 = \frac{30 \, \mu}{\pi \, \delta^4} \; ,
# ```
# ```math
# C_2 = 0 \; ,
# ```
# ```math
# C_3 = \frac{32}{\pi^4 \, \delta^{12}} (\lambda-\mu) \; ,
# ```
# with the first and second Lamé parameters $\lambda$ and $\mu$ and the horizon $\delta$.
#md # !!! danger "Manual definition of the interaction parameters"
#md #     Users can manually specify the interaction parameters $C_1$, $C_2$, and $C_3$.
#md #     The values of the arguments `E` and `nu` are then ignored for the simulation!
#md #     However, for the calculation of `Gc` or `epsilon_c`, they still have to be
#md #     specified when creating a `CPDMaterial` instance!
#md #     Be careful if you use this!

#-
mat = CPDMaterial(horizon=1, rho=8e-6, E=2.1e5, nu=0.25, Gc=2)

#-
mat = CPDMaterial(horizon=1, rho=8e-6, E=2.1e5, nu=0.25, epsilon_c=0.02, C1=1e10)

#src TODO: separate how-to guide!
# ## Multiple material properties within one body
# With the [`MultiMaterial`](@ref) type it is possible to define multiple material
# properties for one single body. Just replace the material definition with a
# `MultiMaterial` instance.
