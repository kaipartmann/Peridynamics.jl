# # [Kalthoff-Winkler experiment](@id tutorial_kalthoff-winkler)
#
# This tutorial demonstrates how to set up and run the Kalthoff-Winkler experiment using the Peridynamics.jl package.
# The Kalthoff-Winkler experiment is a classic dynamic fracture experiment involving a pre-notched sample subjected 
# to impact loading.

# ## Introduction
#
# The Kalthoff-Winkler experiment is widely used to study fracture mechanics under high strain rates.
# This setup provides valuable insights into the behavior of materials under dynamic loading conditions, making it an
# interesting experiment in the field of fracture mechanics.
#
# In this tutorial, we will simulate the Kalthoff-Winkler experiment using peridynamics without the impactor. To do so,
# first we would import the packages.

# ## Import Packages
using Peridynamics

# ## Geometrical Parameters
#
# Define the sample length `l`, width `w`, and thickness `t`, with point spacing `Δx`, horizon size `$ δ $` and crack 
# length `a`.
l  = 200.0E-3  # Length of the sample (meters)
w  = 100.0E-3  # Width of the sample (meters)
t  =   9.0E-3  # Thickness of the sample (meters)
Δx =   1.0E-3  # Discretization size (meters)
δ  = 4.015Δx   # Horizon (meters)
a  = 50.0E-3   # Crack length (meters)

# ## Create the Body
#
# Create a body with the specified dimensions using the bond-based material model with surface corrections:
pos, vol = uniform_box(l, w, t, Δx)
body = Body(BBMaterial{EnergySurfaceCorrection}(), pos, vol)

# ## Material Parameters
#
# The following material parameters are set:
#
# | material parameter | value |
# |:--------|:-------------|
# | Horizon $ δ $ | $4.015 \cdot Δx$ |
# | Young's modulus $E$ | $ 191\cdot 10^{9} \, \mathrm{Pa}$ |
# | Density $ρ$ | $ 8000 ,\mathrm{kg}\,\mathrm{m}^{-3}$ |
# | Griffith's parameter $G_c$ | $22,170 \, \mathrm{J} \, \mathrm{m}^{-2}$ |

material!(body; horizon=4.015Δx, E=191.0e+9, rho=8000.0, Gc=22120.4)

# ## Define Pre-cracks
#
# Define the point sets to insert a crack at the left side of the domain (-x):
point_set!(p -> -a / 2 - δ ≤ p[1] ≤ -a / 2 && 0 ≤ p[2] < w / 2, body, :set_crack1_a)
point_set!(p -> -a / 2 ≤ p[1] ≤ -a / 2 + δ && 0 ≤ p[2] < w / 2, body, :set_crack1_b)
precrack!(body, :set_crack1_a, :set_crack1_b)

# Define the point sets to insert a crack at the right side of the domain (+x):
point_set!(p -> a / 2 - δ ≤ p[1] ≤ a / 2 && 0 ≤ p[2] < w / 2, body, :set_crack2_a)
point_set!(p -> a / 2 ≤ p[1] ≤ a / 2 + δ && 0 ≤ p[2] < w / 2, body, :set_crack2_b)
precrack!(body, :set_crack2_a, :set_crack2_b)

# ## Velocity Boundary Condition
#
# Apply velocity boundary condition on the top edge:
point_set!(p -> -a / 2 + 2.015δ ≤ p[1] ≤ a / 2 - 2.015δ && p[2] ≥ w / 2 - 2.015Δx, body, :set_top)
velocity_bc!(t -> -32.0, body, :set_top, :y)

# ## Simulation
#
# Use the Velocity Verlet algorithm as the time integration method and calculate 2000 time steps:
vv = VelocityVerlet(steps=2000, safety_factor=0.8)

# Define the storage path:
path = joinpath("results", "KW")
ispath(path) && rm(path; recursive=true)  # Delete existing results if they exist

# Create and submit the job:
job = Job(body, vv; path=path, freq=100)
#md # ```julia
#md # @mpitime submit(job)
#md # ```

# ## Visualization
#
# Plot the initial configuration and results (if any).
# Example code to visualize the body and cracks (customize as needed):

# ## Conclusion
#
# This tutorial demonstrated how to set up and run the Kalthoff-Winkler experiment using Peridynamics.jl.
# By simulating this experiment, we can gain insights into the dynamic fracture behavior of materials 
# under high strain rates.
