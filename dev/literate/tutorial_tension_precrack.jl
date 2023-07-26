# # [Mode I tension with predefined crack](@id tutorial_tension_precrack)

# Import the package:
using Peridynamics
using GLMakie #src
using LaTeXStrings #src

#-
# ## Spatial Discretization

# To begin, let's define our geometry. In this example, we will create a point cloud
# consisting of $50 \times 50 \times 5$ material points. The edge length of this block
# is $l = 1\,\text{mm}$, with a thickness of $\frac{1}{10}l$.
#-
l = 1.0
Δx = l / 50
pc = PointCloud(l, l, 0.1l, Δx)

#-
#

#-
fig = Figure(resolution = (1000,500), figure_padding=0) #src
ax = Axis3(fig[1,1]; aspect = :data) #src
hidespines!(ax) #src
hidedecorations!(ax) #src
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:grey80) #src
arrows!(ax, Point3[(0,-l/2,-0.08l),(-l/2,0,-0.08l),(-l/2*1.03,l/2*1.03,0)], #src
        Vec3[(0.1,0,0),(0,0.1,0),(0,0,0.05)], arrowsize=(0.03,0.03,0.05), #src
        linewidth=0.01, color=:black) #src
linesegments!(ax, Point3[(0,-l/2,-0.1l), (0,-l/2,-0.06l), (-l/2,0,-0.1l), #src
              (-l/2,0,-0.06l), (-l/2*1.03-0.02,l/2*1.03,0), #src
              (-l/2*1.03+0.02,l/2*1.03,0)]; linewidth=4, color=:black) #src
text!(ax, L"\mathbf{x}", position=Point3(0.18,-l/2,-0.08l), fontsize=30, #src
      align=(:center,:center), color=:black) #src
text!(ax, L"\mathbf{y}", position=Point3(-l/2,0.18,-0.08l), fontsize=30, #src
      align=(:center,:center), color=:black) #src
text!(ax, L"\mathbf{z}", position=Point3(-l/2*1.03,l/2*1.03,0.13), fontsize=30, #src
      align=(:center,:center), color=:black) #src
save(joinpath(@__DIR__, "..", "assets", "tutorial_tension_1.png"), fig; px_per_unit=3) #src
fig #src

#-
# ![](../assets/tutorial_tension_1.png) #md

#-
# Define a continuum-based material with
# - Horizon $\delta = 3\,\Delta x$
# - Density $\rho = 8e-6\,\mathrm{kg}\,\mathrm{mm}^{-3}$
# - Youngs modulus $E = 210 000 \, \mathrm{MPa}$
# - Poisson ratio $\nu = \frac{1}{4}$
# - Griffith's parameter $G_c = 2.7 \, \mathrm{N} \, \mathrm{mm}^{-1}$
#-
δ = 3.015Δx
mat = ContinuumBasedMaterial(horizon=δ, rho=8e-6, E=2.1e5, nu=0.25, Gc=2.7)

#-
# To add a predefined crack with length $a$, we use two point sets.
#-
a = 0.5l
set_a = findall(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, eachcol(pc.position))
set_b = findall(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, eachcol(pc.position))
precrack = PreCrack(set_a, set_b)

#-
fig = Figure(resolution = (1000,500), figure_padding=0) #src
ax = Axis3(fig[1,1]; aspect = :data) #src
hidespines!(ax) #src
hidedecorations!(ax) #src
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:grey80) #src
arrows!(ax, Point3[(0,-l/2,-0.08l),(-l/2,0,-0.08l),(-l/2*1.03,l/2*1.03,0)], Vec3[(0.1,0,0),(0,0.1,0),(0,0,0.05)], arrowsize=(0.03,0.03,0.05), linewidth=0.01, color=:black) #src
linesegments!(ax, Point3[(0,-l/2,-0.1l), (0,-l/2,-0.06l), (-l/2,0,-0.1l), (-l/2,0,-0.06l), (-l/2*1.03-0.02,l/2*1.03,0), (-l/2*1.03+0.02,l/2*1.03,0)]; linewidth=4, color=:black) #src
text!(ax, L"\mathbf{x}", position=Point3(0.18,-l/2,-0.08l), fontsize=30, align=(:center,:center), color=:black) #src
text!(ax, L"\mathbf{y}", position=Point3(-l/2,0.18,-0.08l), fontsize=30, align=(:center,:center), color=:black) #src
text!(ax, L"\mathbf{z}", position=Point3(-l/2*1.03,l/2*1.03,0.13), fontsize=30, align=(:center,:center), color=:black) #src
meshscatter!(ax, pc.position[:, set_a]; markersize=0.5Δx, color=:blue) #src
meshscatter!(ax, pc.position[:, set_b]; markersize=0.5Δx, color=:red) #src
pca = scatter!(ax, [0,0,0]; markersize=30, color=:blue, label="set_a") #src
pcb = scatter!(ax, [0,0,0]; markersize=30, color=:red, label="set_b") #src
pca.visible[] = false #src
pcb.visible[] = false #src
axislegend(ax; framevisible=false, bgcolor=:transparent, labelsize=26) #src
save(joinpath(@__DIR__, "..", "assets", "tutorial_tension_2.png"), fig; px_per_unit=3) #src
fig #src
# ![](../assets/tutorial_tension_2.png) #md

#-
# Additional point sets for the bottom and the top are used for the velocity boundary condition
# of $\pm 20 \, \mathrm{mm} \, \mathrm{s}^{-1}$.
#-
set_top = findall(p -> p[2] > l/2-Δx, eachcol(pc.position))
set_bottom = findall(p -> p[2] < -l/2+Δx, eachcol(pc.position));

#-
bc_bottom = VelocityBC(t -> -20, set_bottom, 2)

#-
bc_top = VelocityBC(t -> 20, set_top, 2)

#-

#-
fig = Figure(resolution = (1000,500), backgroundcolor = :transparent, figure_padding=0) #src
ax = Axis3(fig[1,1]; aspect = :data) #src
hidespines!(ax) #src
hidedecorations!(ax) #src
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:grey80, transparency=false) #src
arrows!(ax, Point3[(0,-l/2,-0.08l),(-l/2,0,-0.08l),(-l/2*1.03,l/2*1.03,0)], Vec3[(0.1,0,0),(0,0.1,0),(0,0,0.05)], arrowsize=(0.03,0.03,0.05), linewidth=0.01, color=:black) #src
linesegments!(ax, Point3[(0,-l/2,-0.1l), (0,-l/2,-0.06l), (-l/2,0,-0.1l), (-l/2,0,-0.06l), (-l/2*1.03-0.02,l/2*1.03,0), (-l/2*1.03+0.02,l/2*1.03,0)]; linewidth=4, color=:black) #src
text!(ax, L"\mathbf{x}", position=Point3(0.18,-l/2,-0.08l), fontsize=30, align=(:center,:center), color=:black) #src
text!(ax, L"\mathbf{y}", position=Point3(-l/2,0.18,-0.08l), fontsize=30, align=(:center,:center), color=:black) #src
text!(ax, L"\mathbf{z}", position=Point3(-l/2*1.03,l/2*1.03,0.13), fontsize=30, align=(:center,:center), color=:black) #src
meshscatter!(ax, pc.position[:, set_top]; markersize=0.5Δx, color=:blue) #src
meshscatter!(ax, pc.position[:, set_bottom]; markersize=0.5Δx, color=:red) #src
pca = scatter!(ax, [0,0,0]; markersize=30, color=:blue, label="set_top") #src
pcb = scatter!(ax, [0,0,0]; markersize=30, color=:red, label="set_bottom") #src
pca.visible[] = false #src
pcb.visible[] = false #src
axislegend(ax; framevisible=false, bgcolor=:transparent, labelsize=26) #src
save(joinpath(@__DIR__, "..", "assets", "tutorial_tension_3.png"), fig; px_per_unit=3) #src
fig #src
# ![](../assets/tutorial_tension_3.png) #md


#-
# We set the number of time steps for the Velocity Verlet algorithm to 2000 time steps.
vv = VelocityVerlet(2000)

#-
# The results of our analysis should be saved in the directory `"results/mode_I_tension_precrack"` every 10'th time step.
jobname = "mode_I_tension_precrack"
path = joinpath("results", jobname)
!ispath(path) && mkpath(path) # create the path if it does not exist
es = ExportSettings(path, 10)

#-
# Run a single body analysis:
job = PDSingleBodyAnalysis(name=jobname, pc=pc, mat=mat, precracks=[precrack],
                           bcs=[bc_bottom,bc_top], td=vv, es=es)

#-
#md # ```julia
#md # submit(job)
#md # ```

# ### Displacement results:
# ![](../assets/tension_precrack_displacement.gif)

# ### Damage results:
# ![](../assets/tension_precrack_damage.gif)

# (Visualizations made with [ParaView](https://www.paraview.org))
