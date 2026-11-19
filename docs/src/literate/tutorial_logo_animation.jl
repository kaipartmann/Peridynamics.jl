# # [The Peridynamics.jl logo](@id tutorial_logo_animation)

# ```@raw html
#     <video controls loop="true">
#         <source src="https://github.com/user-attachments/assets/4dae9c05-ed70-4ffa-ab4a-b33a74f54531" />
#     </video>
# ```

# (Visualization made with [Makie.jl](https://docs.makie.org))

# The Julia logo crashing into a plate and breaking it into many pieces.
# The last time step of this simulation is the Peridynamics.jl logo.

# First, we have to load the `Peridynamics.jl` package.
using Peridynamics

# ##### Plate
# Now we create the round plate in the background by specifying the diameter, the
# thickness and the point spacing.
lxy = 0.1
lz = 0.005
ΔX₀ₚ = lxy / 100
posₚ, volₚ = round_cylinder(lxy, lz, ΔX₀ₚ)

# To reproduce the logo exactly, the point volumes of the original logo simulation are
# used. There, the volume of a sphere with the diameter of the plate was distributed
# equally on all points of the plate.
volₚ .= 4/3 * π * (lxy / 2)^3 / length(volₚ)
plate = Body(BBMaterial{EnergySurfaceCorrection}(), posₚ, volₚ)

# Then we define the material properties for the plate.
# - Horizon $\delta = 3.015 \Delta x_p$
# - Density $\rho = 2000\,\mathrm{kg}\,\mathrm{m}^{-3}$
# - Youngs modulus $E = 20 \times 10^9 \, \mathrm{Pa}$
# - Griffith's parameter $G_c = 1 \, \mathrm{N} \, \mathrm{m}^{-1}$
material!(plate; horizon=3.015ΔX₀ₚ, E=20e9, rho=2000, Gc=1)

# ##### Julia-logo spheres
# Three spheres are created at the positions of the Julia logo, slightly above the plate.
Ø = 0.03
ΔX₀ₛ = Ø / 40
cz = Ø / 2 + lz / 2 + 1.1 * ΔX₀ₛ
r_logo = Ø / 2 + 0.2 * Ø
sxy, cxy = r_logo * sin(30π / 180), r_logo * cos(30π / 180)
posₛ₁, volₛ₁ = round_sphere(Ø, ΔX₀ₛ; center=(0, r_logo, cz))
posₛ₂, volₛ₂ = round_sphere(Ø, ΔX₀ₛ; center=(cxy, -sxy, cz))
posₛ₃, volₛ₃ = round_sphere(Ø, ΔX₀ₛ; center=(-cxy, -sxy, cz))
sphere₁ = Body(BBMaterial(), posₛ₁, volₛ₁)
sphere₂ = Body(BBMaterial(), posₛ₂, volₛ₂)
sphere₃ = Body(BBMaterial(), posₛ₃, volₛ₃)

# Material properties for the spheres are specified.
# - Horizon $\delta = 3.015 \Delta x_s$
# - Density $\rho = 7850\,\mathrm{kg}\,\mathrm{m}^{-3}$
# - Youngs modulus $E = 210 \times 10^9 \, \mathrm{Pa}$
# - Griffith's parameter $G_c = 10^5 \, \mathrm{N} \, \mathrm{m}^{-1}$

# All material points of the spheres have a initial velocity
# of $-30\, \mathrm{m} \, \mathrm{s}^{-1}$ in $z$-direction.
# The spheres should not break, so failure is deactivated for them.
for sphere in (sphere₁, sphere₂, sphere₃)
    material!(sphere; horizon=3.015ΔX₀ₛ, E=210e9, rho=7850, Gc=1e5)
    velocity_ic!(sphere, :all_points, :z, -30)
    no_failure!(sphere)
end

# For the contact analysis, all bodies need to be specified in a [`MultibodySetup`](@ref).
ms = MultibodySetup(:plate => plate, :sphere1 => sphere₁, :sphere2 => sphere₂,
                    :sphere3 => sphere₃)

# Contact between the plate and the three spheres needs to be specified.
contact!(ms, :plate, :sphere1; radius=1.2ΔX₀ₚ, penalty_factor=1.2e12)
contact!(ms, :plate, :sphere2; radius=1.2ΔX₀ₚ, penalty_factor=1.2e12)
contact!(ms, :plate, :sphere3; radius=1.2ΔX₀ₚ, penalty_factor=1.2e12)

# The simulation time is $0.6\,\mathrm{ms}$ with explicit time integration using the
# Velocity Verlet algorithm, which results in 10250 time steps.
vv = VelocityVerlet(time=0.0006)

# Now we create a directory for the results and create a Job.
# The results are exported every 50 time steps.
job = Job(ms, vv; path="results/logo", freq=50)

# To complete everything, the Job is submitted for simulation.
#md # ```julia
#md # submit(job)
#md # ```

# ##### Animation
# To create the animation, every exported time step is plotted with
# [GLMakie.jl](https://docs.makie.org) and saved as an image.
#md # ```julia
#md # using GLMakie
#md # ```

# The following function plots one time step. The plate is colored by its damage, the
# spheres get the colors of the Julia logo. Points of the plate that are completely broken
# off (damage > 96%) are hidden. The size of the figure and the camera are set to the same
# view as in the logo.
#md # ```julia
#md # function plot_logo(step)
#md #     vtk(body) = read_vtk("results/logo/vtk/$(body)_timestep_$(step).pvtu")
#md #     plate = vtk("plate")
#md #     ids = findall(d -> d < 0.96, plate[:damage])
#md #     colors = Makie.Colors.JULIA_LOGO_COLORS
#md #     damage_colors = [RGBf(0.251, 0.338, 0.847), RGBf(0.865, 0.865, 0.865),
#md #                      RGBf(0.706, 0.016, 0.149)]
#md #
#md #     fig = Figure(size=(674, 678), figure_padding=0)
#md #     ax = LScene(fig[1, 1]; show_axis=false)
#md #     cam3d!(ax.scene; center=false, eyeposition=Vec3d(-0.0207, -0.00288, 0.142),
#md #            lookat=Vec3d(-0.00218, -0.00025, 0), upvector=Vec3d(0.0656, 0.9975, 0.027),
#md #            fov=40.62)
#md #     meshscatter!(ax, plate[:position][:, ids]; markersize=0.6ΔX₀ₚ,
#md #                  color=plate[:damage][ids], colormap=damage_colors, colorrange=(0, 1))
#md #     meshscatter!(ax, vtk("sphere1")[:position]; markersize=0.6ΔX₀ₛ,
#md #                  color=colors.green)
#md #     meshscatter!(ax, vtk("sphere2")[:position]; markersize=0.6ΔX₀ₛ,
#md #                  color=colors.purple)
#md #     meshscatter!(ax, vtk("sphere3")[:position]; markersize=0.6ΔX₀ₛ,
#md #                  color=colors.red)
#md #     return fig
#md # end
#md # ```

# Now every exported time step is plotted and saved as an image.
# The last image is the Peridynamics.jl logo.
#md # ```julia
#md # mkpath("results/logo/frames")
#md # for (i, step) in enumerate(0:50:10250)
#md #     save("results/logo/frames/frame_$(lpad(i, 4, '0')).png", plot_logo(step);
#md #          px_per_unit=2)
#md # end
#md # ```

# Finally, the images are combined into a video with [FFMPEG](https://ffmpeg.org),
# which is installed together with Makie.jl.
#md # ```julia
#md # ffmpeg = Makie.FFMPEG_jll.ffmpeg()
#md # run(`$ffmpeg -y -loglevel error -framerate 25 -i results/logo/frames/frame_%04d.png
#md #      -pix_fmt yuv420p results/logo/logo.mp4`)
#md # ```

# ##### Visualization with ParaView
# Alternatively, the results can be visualized with [ParaView](https://www.paraview.org),
# as explained in [Visualization with ParaView](@ref visualization).
# To get the same view as in the logo:
# 1. Open all `.pvtu` files in `results/logo/vtk` and click **Apply**.
# 2. Hide the broken points of the plate with a **Threshold** filter on `damage` from 0 to
#    0.96.
# 3. Show the points as spheres: representation **Point Gaussian**, shader preset
#    *Sphere*, and a **Gaussian Radius** of `0.0006` for the plate and `0.00045` for the
#    spheres.
# 4. Color the plate by `damage` and the spheres green, purple and red.
# 5. Set the camera with **Adjust Camera**:
#    - Position: `-0.0207, -0.00288, 0.142`
#    - Focal Point: `-0.00218, -0.00025, 0`
#    - View Up: `0.0656, 0.9975, 0.027`
#    - View Angle: `40.62`
# 6. Save the images or a video with **File** $\rightarrow$ **Save Animation**.
