# # [Tension with predefined crack](@id tutorial_tension_precrack)

# Import the package:
using Peridynamics
# using GLMakie #src
# using LaTeXStrings #src

# First some geometrical parameters are defined.
# These are edge length `l`, point spacing `Δx` and crack length `a`.
# Now a cuboid body with the specified edge lengths and a thickness of one tenth thereof is
# created using the bond-based material model.
l, Δx, a = 1.0, 1/50, 0.5
pos, vol = uniform_box(l, l, 0.1l, Δx)
body = Body(BBMaterial(), pos, vol)

#-
# fig = Figure(resolution = (1000,500), figure_padding=0) #src
# ax = Axis3(fig[1,1]; aspect = :data) #src
# hidespines!(ax) #src
# hidedecorations!(ax) #src
# meshscatter!(ax, pc.position; markersize=0.5Δx, color=:grey80) #src
# arrows!(ax, Point3[(0,-l/2,-0.08l),(-l/2,0,-0.08l),(-l/2*1.03,l/2*1.03,0)], #src
#         Vec3[(0.1,0,0),(0,0.1,0),(0,0,0.05)], arrowsize=(0.03,0.03,0.05), #src
#         linewidth=0.01, color=:black) #src
# linesegments!(ax, Point3[(0,-l/2,-0.1l), (0,-l/2,-0.06l), (-l/2,0,-0.1l), #src
#               (-l/2,0,-0.06l), (-l/2*1.03-0.02,l/2*1.03,0), #src
#               (-l/2*1.03+0.02,l/2*1.03,0)]; linewidth=4, color=:black) #src
# text!(ax, L"\mathbf{x}", position=Point3(0.18,-l/2,-0.08l), fontsize=30, #src
#       align=(:center,:center), color=:black) #src
# text!(ax, L"\mathbf{y}", position=Point3(-l/2,0.18,-0.08l), fontsize=30, #src
#       align=(:center,:center), color=:black) #src
# text!(ax, L"\mathbf{z}", position=Point3(-l/2*1.03,l/2*1.03,0.13), fontsize=30, #src
#       align=(:center,:center), color=:black) #src
# # save(joinpath(@__DIR__, "..", "assets", "tutorial_tension_1.png"), fig; px_per_unit=3) #src
# fig #src

#-
# ![](https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/ca3e724b-9dd7-417a-b3aa-bc37f81f4674) #md

# The following material parameters are set:
#
# | material parameter | value |
# |:--------|:-------------|
# | Horizon $ δ $ | $3.015 \cdot Δx$ |
# | Young's modulus $E$ | $ 210000 \, \mathrm{MPa}$ |
# | Density $ρ$ | $ 8 \cdot 10^{-6}\,\mathrm{kg}\,\mathrm{mm}^{-3}$ |
# | Griffith's parameter $G_c$ | $2.7 \, \mathrm{N} \, \mathrm{mm}^{-1}$ |
δ = 3.015Δx
material!(body; horizon=δ, E=2.1e5, rho=8e-6, Gc=2.7)

# Two point sets are defined to insert a crack between them:
point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, body, :set_a)
point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, body, :set_b)
precrack!(body, :set_a, :set_b)

# fig = Figure(resolution = (1000,500), figure_padding=0) #src
# ax = Axis3(fig[1,1]; aspect = :data) #src
# hidespines!(ax) #src
# hidedecorations!(ax) #src
# meshscatter!(ax, pc.position; markersize=0.5Δx, color=:grey80) #src
# arrows!(ax, Point3[(0,-l/2,-0.08l),(-l/2,0,-0.08l),(-l/2*1.03,l/2*1.03,0)], Vec3[(0.1,0,0),(0,0.1,0),(0,0,0.05)], arrowsize=(0.03,0.03,0.05), linewidth=0.01, color=:black) #src
# linesegments!(ax, Point3[(0,-l/2,-0.1l), (0,-l/2,-0.06l), (-l/2,0,-0.1l), (-l/2,0,-0.06l), (-l/2*1.03-0.02,l/2*1.03,0), (-l/2*1.03+0.02,l/2*1.03,0)]; linewidth=4, color=:black) #src
# text!(ax, L"\mathbf{x}", position=Point3(0.18,-l/2,-0.08l), fontsize=30, align=(:center,:center), color=:black) #src
# text!(ax, L"\mathbf{y}", position=Point3(-l/2,0.18,-0.08l), fontsize=30, align=(:center,:center), color=:black) #src
# text!(ax, L"\mathbf{z}", position=Point3(-l/2*1.03,l/2*1.03,0.13), fontsize=30, align=(:center,:center), color=:black) #src
# meshscatter!(ax, pc.position[:, set_a]; markersize=0.5Δx, color=:blue) #src
# meshscatter!(ax, pc.position[:, set_b]; markersize=0.5Δx, color=:red) #src
# pca = scatter!(ax, [0,0,0]; markersize=30, color=:blue, label="set_a") #src
# pcb = scatter!(ax, [0,0,0]; markersize=30, color=:red, label="set_b") #src
# pca.visible[] = false #src
# pcb.visible[] = false #src
# axislegend(ax; framevisible=false, bgcolor=:transparent, labelsize=26) #src
# save(joinpath(@__DIR__, "..", "assets", "tutorial_tension_2.png"), fig; px_per_unit=3) #src
# fig #src
# ![](https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/36a16d3e-a2a5-4db7-9b73-55b321583902) #md

# Two more point sets at the top and at the bottom are created, which are used for the
# velocity boundary condition.
point_set!(p -> p[2] > l/2-Δx, body, :set_top)
point_set!(p -> p[2] < -l/2+Δx, body, :set_bottom)

# fig = Figure(resolution = (1000,500), backgroundcolor = :transparent, figure_padding=0) #src
# ax = Axis3(fig[1,1]; aspect = :data) #src
# hidespines!(ax) #src
# hidedecorations!(ax) #src
# meshscatter!(ax, pc.position; markersize=0.5Δx, color=:grey80, transparency=false) #src
# arrows!(ax, Point3[(0,-l/2,-0.08l),(-l/2,0,-0.08l),(-l/2*1.03,l/2*1.03,0)], Vec3[(0.1,0,0),(0,0.1,0),(0,0,0.05)], arrowsize=(0.03,0.03,0.05), linewidth=0.01, color=:black) #src
# linesegments!(ax, Point3[(0,-l/2,-0.1l), (0,-l/2,-0.06l), (-l/2,0,-0.1l), (-l/2,0,-0.06l), (-l/2*1.03-0.02,l/2*1.03,0), (-l/2*1.03+0.02,l/2*1.03,0)]; linewidth=4, color=:black) #src
# text!(ax, L"\mathbf{x}", position=Point3(0.18,-l/2,-0.08l), fontsize=30, align=(:center,:center), color=:black) #src
# text!(ax, L"\mathbf{y}", position=Point3(-l/2,0.18,-0.08l), fontsize=30, align=(:center,:center), color=:black) #src
# text!(ax, L"\mathbf{z}", position=Point3(-l/2*1.03,l/2*1.03,0.13), fontsize=30, align=(:center,:center), color=:black) #src
# meshscatter!(ax, pc.position[:, set_top]; markersize=0.5Δx, color=:blue) #src
# meshscatter!(ax, pc.position[:, set_bottom]; markersize=0.5Δx, color=:red) #src
# pca = scatter!(ax, [0,0,0]; markersize=30, color=:blue, label="set_top") #src
# pcb = scatter!(ax, [0,0,0]; markersize=30, color=:red, label="set_bottom") #src
# pca.visible[] = false #src
# pcb.visible[] = false #src
# axislegend(ax; framevisible=false, bgcolor=:transparent, labelsize=26) #src
# save(joinpath(@__DIR__, "..", "assets", "tutorial_tension_3.png"), fig; px_per_unit=3) #src
# fig #src
# ![](https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/7c3e0b71-57b7-4bff-9bca-617bdae9f5be) #md

# The tension is applied by moving the ends of the body apart at a constant speed of
# $\pm 50 \, \mathrm{mm} \, \mathrm{s}^{-1}$:
velocity_bc!(t -> -30, body, :set_bottom, :y)
velocity_bc!(t -> 30, body, :set_top, :y)
# The Velocity Verlet algorithm is used as time integration method and 2000 time steps
# are calculated:
vv = VelocityVerlet(steps=2000)
# Now the job is defined
job = Job(body, vv; path="results/mode_i_tension_precrack")
# Finally the job is submitted to start simulations
#md # ```julia
#md # submit(job)
#md # ```

# ### Damage results:

# ```@raw html
#     <video controls loop="true">
#         <source src="https://github.com/user-attachments/assets/918aed7b-735d-418f-900c-e0a996db2bab" />
#     </video>
# ```
