using Peridynamics

include("sphere_shape_trixip.jl")

function new_logo(path; np=50, ns=15)
    # plate
    lxy = 0.1
    lz = 0.005
    ΔX₀ₚ = lxy / np
    posₚ, volₚ = cylinder(lxy, lz, ΔX₀ₚ)
    plate = Body(BBMaterial{EnergySurfaceCorrection}(), posₚ, volₚ)
    material!(plate; horizon=3.015ΔX₀ₚ, E=20e9, rho=2000, Gc=1)

    # spheres
    Ø = 0.03
    r = Ø / 2
    ΔX₀ₛ = Ø / ns
    cz = Ø / 2 + lz / 2 + 1.1 * ΔX₀ₛ
    r_logo = Ø / 2 + 0.2 * Ø
    sxy, cxy = r_logo * sin(30π / 180), r_logo * cos(30π / 180)
    posₛ₁, volₛ₁ = sphere_shape(ΔX₀ₛ, r, (0.0, r_logo, cz))
    posₛ₂, volₛ₂ = sphere_shape(ΔX₀ₛ, r, (cxy, -sxy, cz))
    posₛ₃, volₛ₃ = sphere_shape(ΔX₀ₛ, r, (-cxy, -sxy, cz))
    sphere₁ = Body(BBMaterial(), posₛ₁, volₛ₁)
    sphere₂ = Body(BBMaterial(), posₛ₂, volₛ₂)
    sphere₃ = Body(BBMaterial(), posₛ₃, volₛ₃)
    ms = MultibodySetup(:plate => plate, :sphere1 => sphere₁, :sphere2 => sphere₂,
                        :sphere3 => sphere₃)
    for sphere in (sphere₁, sphere₂, sphere₃)
        material!(sphere; horizon=3.015ΔX₀ₛ, E=210e9, rho=7850, Gc=1e5)
        velocity_ic!(sphere, :all_points, :z, -30)
        failure_permit!(sphere, false)
    end
    contact!(ms, :plate, :sphere1; radius=1.2ΔX₀ₚ, penalty_factor=1.2e12)
    contact!(ms, :plate, :sphere2; radius=1.2ΔX₀ₚ, penalty_factor=1.2e12)
    contact!(ms, :plate, :sphere3; radius=1.2ΔX₀ₚ, penalty_factor=1.2e12)

    vv = VelocityVerlet(time=0.0006)
    job = Job(ms, vv; path=path)
    return job
end

##
path = "results/newlogo_v2"
# rm(path; recursive=true, force=true)
np, ns = 30, 6
job = new_logo(path, np=np, ns=ns)
@time submit(job)
np, ns = 100, 40
job = new_logo(path, np=np, ns=ns)
@time submit(job)

##
using CairoMakie

rad(grad::Real) = grad * π / 180
# paraview = cgrad([RGBf(0.231373, 0.298039, 0.752941),
#                   RGBf(0.865003, 0.865003, 0.865003),
#                   RGBf(0.705882, 0.0156863, 0.14902)])
dmg = cgrad([RGBf(0.251, 0.338, 0.847),
             RGBf(0.865003, 0.865003, 0.865003),
             RGBf(0.705882, 0.0156863, 0.14902)])

##
r_p = 0.1 / np * 0.6
r_s = 0.03 / ns * 0.6

##
timestep = "010250"
sim = "newlogo_v2"
vtk_file_p = "results/$sim/vtk/plate_timestep_$timestep.pvtu"
vtk_file_s1 = "results/$sim/vtk/sphere1_timestep_$timestep.pvtu"
vtk_file_s2 = "results/$sim/vtk/sphere2_timestep_$timestep.pvtu"
vtk_file_s3 = "results/$sim/vtk/sphere3_timestep_$timestep.pvtu"
res_p = read_vtk(vtk_file_p)
res_s1 = read_vtk(vtk_file_s1)
res_s2 = read_vtk(vtk_file_s2)
res_s3 = read_vtk(vtk_file_s3)

pids = findall(x -> x < 0.96, res_p[:damage])

fig = Figure(size=(1000, 1000), figure_padding=5)
ax = LScene(fig[1, 1]; show_axis=false)
cam3d_cad!(ax.scene)
cam3d_cad!(fig.scene)
ms_p = meshscatter!(ax, res_p[:position][:,pids]; markersize=r_p,
                    color=res_p[:damage][pids], colormap=dmg, colorrange=(0, 1))
ms_s1 = meshscatter!(ax, res_s1[:position]; markersize=r_s, color=RGBf(0.22, 0.596, 0.149))
ms_s2 = meshscatter!(ax, res_s2[:position]; markersize=r_s, color=RGBf(0.584, 0.345, 0.693))
ms_s3 = meshscatter!(ax, res_s3[:position]; markersize=r_s, color=RGBf(0.796, 0.235, 0.2))
# Colorbar(fig[1,2], ms_p; label="damage", height=120, tellheight=false)
colgap!(fig.layout, 0)
# rotate_cam!(ax.scene, 0, rad(210), 0)
# rotate_cam!(ax.scene, rad(-50), rad(-1), rad(10))
# rotate_cam!(ax.scene, rad(-30), 0, 0)
# rotate_cam!(ax.scene, 0, 0, 0)

cam = cameracontrols(ax.scene)
cam.lookat[] = [0, 0, 0]
cam.eyeposition[] = [-0.1, -0.01, 0.5]
update_cam!(ax.scene, cam) # not sure if there is a method without cam

rotate_cam!(ax.scene, 0, 0, rad(-10))
fig
