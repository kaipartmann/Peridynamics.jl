# # [Kalthoff-Winkler experiment](@id tutorial_kalthoff-winkler)

# Import the package:
using Peridynamics
using GLMakie #src
using LaTeXStrings #src

# This tutorial demonstrates how to set up and run the Kalthoff-Winkler experiment using the Peridynamics.jl package.
# The Kalthoff-Winkler experiment is a classic dynamic fracture experiment involving a pre-notched sample subjected to an impact loading.

# The geometrical parameters are defined as:
# Sample length `l`, width `w`, and thickness `t`, with point spacing `Δx`, horizon size `$ δ $` and crack length `a`.
l, w, t, Δx, δ, a = 100.0, 100 / 50, 50.0

# Now a cuboid body with the specified dimensions is created using the bond-based material model with surface corrections:
pos, vol = uniform_box(l, w, t, Δx)
body = Body(BBMaterial{EnergySurfaceCorrection}(), pos, vol)

# Some plotting

# The following material parameters are set:
#
# | material parameter | value |
# |:--------|:-------------|
# | Horizon $ δ $ | $4.015 \cdot Δx$ |
# | Young's modulus $E$ | $ 191\cdot 10^{9} \, \mathrm{Pa}$ |
# | Density $ρ$ | $ 8000 ,\mathrm{kg}\,\mathrm{m}^{-3}$ |
# | Griffith's parameter $G_c$ | $22120.4 \, \mathrm{N} \, \mathrm{m}^{-1}$ |

material!(body; horizon=4.015Δx, E=191.0e+9, rho=8000.0, Gc=22120.4)

# Two point sets are defined to insert a crack between them at the left side of the domain (-x):
point_set!(p -> -a / 2 - δ ≤ p[1] ≤ -a / 2 && 0 ≤ p[2] < wid / 2, body, :set_crack1_a)
point_set!(p -> -a / 2 ≤ p[1] ≤ -a / 2 + δ && 0 ≤ p[2] < wid / 2, body, :set_crack1_b)
precrack!(body, :set_crack1_a, :set_crack1_b)

# For the second crack, the same points are defined but on the right side of the domain (+x):
point_set!(p -> a / 2 - δ ≤ p[1] ≤ a / 2 && 0 ≤ p[2] < wid / 2, body, :set_crack2_a)
point_set!(p -> a / 2 ≤ p[1] ≤ a / 2 + δ && 0 ≤ p[2] < wid / 2, body, :set_crack2_b)
precrack!(body, :set_crack2_a, :set_crack2_b)

# A point sets at the top is created, which is used for the velocity boundary condition:
point_set!(p -> -a / 2 + 2.015δ ≤ p[1] ≤ a / 2 - 2.015δ && p[2] ≥ wid / 2 - 2.015Δx, body,
           :set_top)


# The specimen is impacted by a projectile moving at a speed of $\pm 32 \, \mathrm{m} \, \mathrm{s}^{-1}$:
velocity_bc!(t -> -32.0, body, :set_top, :y)

# The Velocity Verlet algorithm is used as time integration method and 2000 time steps
# are calculated:
vv = VelocityVerlet(steps=2000, safety_factor=0.8)
# Then the storage path is defined:
root = joinpath(@__DIR__, "results", "mode_i")
#  Now the job is defined
job = Job(body, vv; path=path, freq=100)
# Finally the job is submitted to start simulations
#md # ```julia
#md # submit(job)
#md # ```
