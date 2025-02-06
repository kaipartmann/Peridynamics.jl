# # Wave propagation across material interface

# Based on the [wave propagation tutorial](@ref tut_wave), this tutorial features a
# displacement wave crossing a material interface, as investigated in
# [Partmann2024AAM](@cite). This means that a bar containing two sections with different
# material properties is regarded.

# First import the Peridynamics.jl package:
using Peridynamics

# Then the geometric parameters are set, which are the same as in the wave propagation
# tutorial.
# These are the length `lx`, width and height `lyz` of the cuboid as well as the number
# of points in the width `npyz`, which determines the point spacing `Δx`.
lx = 0.2
lyz = 0.002
npyz = 4
Δx = lyz / npyz
# With these parameters we now create a body, here using the non-ordinary state-based
# correspondence formulation.
pos, vol = uniform_box(lx, lyz, lyz, Δx)
body = Body(NOSBMaterial(), pos, vol)
# Again, failure is not allowed in the whole body.
no_failure!(body)

# Then the material parameters for one half of the body are assigned to the whole body
# first.
material!(body, horizon=3.015Δx, rho=7850.0, E=210e9, nu=0.25, epsilon_c=0.01)

# Now a point set containing the other half of all points is created.
point_set!(x -> x < 0, body, :set1)
# ![](https://github.com/user-attachments/assets/c589ab03-e7f5-4b71-b0c2-0f05eb3cdddc)
# The parameters for this point set are then overwritten with their new parameters.
material!(body, :set1, horizon=3.015Δx, rho=7850.0, E=105e9, nu=0.25, epsilon_c=0.01)

# Except for the Young's modulus, these are the same in both sections:
#
# | material parameter | value |
# |:--------|:-------------|
# | Horizon $ δ $ | $3.015 \cdot Δx$ |
# | Density $ρ$ | $ 7850\,\mathrm{kg}\,\mathrm{m}^{-3}$ |
# | Young's modulus $E_I$ | $ 105 \, \mathrm{GPa}$ |
# | Young's modulus $E_{II}$ | $ 210 \, \mathrm{GPa}$ |
# | Poisson's ratio $ν$ | $0.25$ |
# | critical strain $ε_c$ | $0.01$ |


# To employ the boundary conditions creating a displacement wave, the point set `:left` is
# created:
point_set!(x -> x < -lx / 2 + 1.2Δx, body, :left)
# As in the wave propagation tutorial, the applied velocity boundary condition is
# ```math
#     {v}_x (t) =
#     \begin{cases}
#         v_\mathrm{max} \cdot \sin(2\pi \cdot \frac{t}{T}) \qquad
#         & \forall \; 0 \leq t \leq T \\
#         0 &\text{else.}
#     \end{cases}
# ```
# ![](https://github.com/user-attachments/assets/34165b09-2df9-495d-8e41-acff28d0098c)
T, vmax = 1.0e-5, 2.0
velocity_bc!(t -> t < T ? vmax * sin(2π / T * t) : 0.0, body, :left, :x)

# The Velocity Verlet algorithm is used as time integration method and 2000 time steps
# are calculated:
vv = VelocityVerlet(steps=2000)
# Finally the job is defined and submitted.
#md # ```julia
#md # job = Job(body, vv; path="results/xwave_interface")
#md # submit(job)
#md # ```

# The following video shows the results for the displacement in x-direction $u_x$ and the
# first entry of the
# stress tensor, i.e. the normal stress in x-direction $\sigma_{11}$.
# ```@raw html
#     <video controls loop="true">
#         <source src="https://github.com/user-attachments/assets/d7e600ee-cad8-40f9-99f4-e96ac05585d6" />
#     </video>
# ```
