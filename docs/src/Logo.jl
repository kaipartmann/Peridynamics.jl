using Peridynamics

## Setup
# plate
lxy = 0.1
lz = 0.005
ΔX₀ₚ = lxy / 100
posₚ, volₚ = round_cylinder(lxy, lz, ΔX₀ₚ)
plate = Body(BBMaterial{EnergySurfaceCorrection}(), posₚ, volₚ)
material!(plate; horizon=3.015ΔX₀ₚ, E=20e9, rho=2000, Gc=1)

# spheres
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
for sphere in (sphere₁, sphere₂, sphere₃)
    material!(sphere; horizon=3.015ΔX₀ₛ, E=210e9, rho=7850, Gc=1e5)
    velocity_ic!(sphere, :all_points, :z, -30)
    no_failure!(sphere)
end

ms = MultibodySetup(:plate => plate, :sphere1 => sphere₁, :sphere2 => sphere₂,
                    :sphere3 => sphere₃)
contact!(ms, :plate, :sphere1; radius=1.2ΔX₀ₚ, penalty_factor=1.2e12)
contact!(ms, :plate, :sphere2; radius=1.2ΔX₀ₚ, penalty_factor=1.2e12)
contact!(ms, :plate, :sphere3; radius=1.2ΔX₀ₚ, penalty_factor=1.2e12)

# the last time step is the frame of the logo
vv = VelocityVerlet(steps=2000)
job = Job(ms, vv; path="results/logo", freq=10)

## Simulation
@time submit(job)

## Plot every exported time step
using GLMakie

function plot_logo(step)
    vtk(body) = read_vtk("results/logo/vtk/$(body)_timestep_$(step).pvtu")
    plate = vtk("plate")
    ids = findall(d -> d < 0.96, plate[:damage])
    colors = Makie.Colors.JULIA_LOGO_COLORS
    damage_colors = [RGBf(0.251, 0.338, 0.847), RGBf(0.865, 0.865, 0.865),
                     RGBf(0.706, 0.016, 0.149)]

    fig = Figure(size=(1000, 1000), figure_padding=0)
    ax = LScene(fig[1, 1]; show_axis=false)
    cam3d!(ax.scene; center=false, eyeposition=Vec3d(-0.0166, -0.003, 0.0929),
           lookat=Vec3d(-0.0034, -0.0007, 0), upvector=Vec3d(0.07, 0.997, 0.0342),
           fov=58)
    meshscatter!(ax, plate[:position][:, ids]; markersize=0.6ΔX₀ₚ,
                 color=plate[:damage][ids], colormap=damage_colors, colorrange=(0, 1))
    meshscatter!(ax, vtk("sphere1")[:position]; markersize=0.6ΔX₀ₛ,
                 color=colors.green)
    meshscatter!(ax, vtk("sphere2")[:position]; markersize=0.6ΔX₀ₛ,
                 color=colors.purple)
    meshscatter!(ax, vtk("sphere3")[:position]; markersize=0.6ΔX₀ₛ,
                 color=colors.red)
    return fig
end

mkpath("results/logo/frames")
for (i, step) in enumerate(0:10:2000)
    save("results/logo/frames/frame_$(lpad(i, 4, '0')).png", plot_logo(step))
end

## Video
ffmpeg = Makie.FFMPEG_jll.ffmpeg()
run(`$ffmpeg -y -loglevel error -framerate 25 -i results/logo/frames/frame_%04d.png
     -pix_fmt yuv420p results/logo/logo.mp4`)
