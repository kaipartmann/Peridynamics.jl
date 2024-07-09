# # [The old Peridynamics.jl logo](@id tutorial_logo)

# ![](https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/b79fd475-67fc-4097-9f05-791833f29fa9)

# (Visualization made with [ParaView](https://www.paraview.org))

# The Julia logo crashing into a plate and braking it into many pieces.

# First, we have to load the `Peridynamics.jl` package.
using Peridynamics

# ##### Plate
# Now we create the plate in the background by specifying the dimensions and the point
# spacing.
lxy = 0.1
lz = 0.01
ΔX₀ₚ = lxy / 50
posₚ, volₚ = uniform_box(lxy, lxy, lz, ΔX₀ₚ)
plate = Body(BBMaterial{EnergySurfaceCorrection}(), posₚ, volₚ)

# Then we define the material properties for the plate.
# - Horizon $\delta = 3.015 \Delta x_p$
# - Density $\rho = 2000\,\mathrm{kg}\,\mathrm{m}^{-3}$
# - Youngs modulus $E = 30 \times 10^9 \, \mathrm{Pa}$
# - Griffith's parameter $G_c = 10 \, \mathrm{N} \, \mathrm{m}^{-1}$
material!(plate; horizon=3.015ΔX₀ₚ, E=30e9, rho=2000, Gc=10)

# ##### Julia-logo spheres
# A spherical body is created, where only the points inside a specified radius are
# preserved to create the spheres of the logo.
# These points are then copied three times and moved to the correct position to represent
# the logo.
Ø = 0.03
ΔX₀ₛ = Ø / 20
cz = Ø / 2 + lz / 2 + 1.1 * ΔX₀ₛ
r_logo = Ø / 2 + 0.2 * Ø
sxy, cxy = r_logo * sin(30π / 180), r_logo * cos(30π / 180)
posₛ₁, volₛ₁ = uniform_sphere(Ø, ΔX₀ₛ; center_y=r_logo, center_z=cz)
posₛ₂, volₛ₂ = uniform_sphere(Ø, ΔX₀ₛ; center_x=cxy, center_y=-sxy, center_z=cz)
posₛ₃, volₛ₃ = uniform_sphere(Ø, ΔX₀ₛ; center_x=-cxy, center_y=-sxy, center_z=cz)
sphere₁ = Body(BBMaterial(), posₛ₁, volₛ₁)
sphere₂ = Body(BBMaterial(), posₛ₂, volₛ₂)
sphere₃ = Body(BBMaterial(), posₛ₃, volₛ₃)

# Material properties for the spheres are specified.
# - Horizon $\delta = 3.015 \Delta x_s$
# - Density $\rho = 7850\,\mathrm{kg}\,\mathrm{m}^{-3}$
# - Youngs modulus $E = 210 \times 10^9 \, \mathrm{Pa}$
# - Griffith's parameter $G_c = 1000 \, \mathrm{N} \, \mathrm{m}^{-1}$

# All material points of the spheres have a initial velocity
# of $-20\, \mathrm{m} \, \mathrm{s}^{-1}$ in $z$-direction.
for sphere in (sphere₁, sphere₂, sphere₃)
    material!(sphere; horizon=3.015ΔX₀ₛ, E=210e9, rho=7850, Gc=1000)
    velocity_ic!(sphere, :all_points, :z, -20)
end

# Multibody Setup?
# For the contact analysis, all bodies need to be specified in a [`MultibodySetup`](@ref).
ms = MultibodySetup(:plate => plate, :sphere1 => sphere₁, :sphere2 => sphere₂,
                    :sphere3 => sphere₃)

# Contact between the plate and the three spheres needs to be specified.
contact!(ms, :plate, :sphere1; radius=ΔX₀ₚ)
contact!(ms, :plate, :sphere2; radius=ΔX₀ₚ)
contact!(ms, :plate, :sphere3; radius=ΔX₀ₚ)

# For this simulation, 3000 time steps with explicit time integration and the velocity
# verlet algorithm are used.
vv = VelocityVerlet(steps=3000)

# Now we create a directory for the results and create a Job.
job = Job(ms, vv; path="results/logo")

rm(path; recursive=true, force=true) #src

# To complete everything, the Job is submitted for simulation.
#md # ```julia
#md # submit(job)
#md # ```
