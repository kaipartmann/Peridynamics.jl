using Peridynamics #hide
using GLMakie #src

#-
# # [Predefined cracks](@id precracks)
# Initially defined cracks can be defined with [`PreCrack`](@ref).
# First, we define a point cloud:
lx = 10
ly = 10
lz = 1
Δx = 0.5
pc = PointCloud(lx, ly, lz, Δx)

#-
fig = Figure(resolution = (1000,700), backgroundcolor = :transparent) #src
ax = Axis3(fig[1,1]; aspect = :data, elevation = 0.4π, azimuth = 1.45π) #src
enlfac = 1.05 #src
xlims!(ax, (-lx/2 * enlfac, lx/2 * enlfac)) #src
ylims!(ax, (-ly/2 * enlfac, ly/2 * enlfac)) #src
zlims!(ax, (-lz/2 * enlfac, lz/2 * enlfac)) #src
hidespines!(ax) #src
hidedecorations!(ax) #src
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:grey70) #src
save(joinpath(@__DIR__, "..", "assets", "pc_precrack.png"), fig; px_per_unit=3) #src
fig #src
# ![](../assets/pc_precrack.png) #md

#-
# Then, we define a crack length and specify the point sets of points that should not
# interact with each other.

#-
cracklength = 0.5 * lx

#-
set_a = findall(
    (pc.position[2, :] .>= 0) .&
    (pc.position[2, :] .< 6Δx) .&
    (pc.position[1, :] .<= -lx/2 + cracklength),
)

#-
set_b = findall(
    (pc.position[2, :] .<= 0) .&
    (pc.position[2, :] .> -6Δx) .&
    (pc.position[1, :] .<= -lx/2 + cracklength),
)

#-
fig = Figure(resolution = (1000,700), backgroundcolor = :transparent) #src
ax = Axis3(fig[1,1]; aspect = :data, elevation = 0.4π, azimuth = 1.45π) #src
enlfac = 1.05 #src
xlims!(ax, (-lx/2 * enlfac, lx/2 * enlfac)) #src
ylims!(ax, (-ly/2 * enlfac, ly/2 * enlfac)) #src
zlims!(ax, (-lz/2 * enlfac, lz/2 * enlfac)) #src
hidespines!(ax) #src
hidedecorations!(ax) #src
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:grey70) #src
meshscatter!(ax, pc.position[:, set_a]; markersize=0.5Δx, color=:blue) #src
meshscatter!(ax, pc.position[:, set_b]; markersize=0.5Δx, color=:red) #src
pca = scatter!(ax, [0,0,0]; markersize=30, color=:blue, label="set_a") #src
pcb = scatter!(ax, [0,0,0]; markersize=30, color=:red, label="set_b") #src
pca.visible[] = false #src
pcb.visible[] = false #src
Legend(fig[1,2], ax; framevisible=false, bgcolor=:transparent, labelsize=28) #src
save(joinpath(@__DIR__, "..", "assets", "pc_precrack_with_sets.png"), fig; px_per_unit=3) #src
fig #src
# ![](../assets/pc_precrack_with_sets.png) #md

#-
precrack = PreCrack(set_a, set_b)
