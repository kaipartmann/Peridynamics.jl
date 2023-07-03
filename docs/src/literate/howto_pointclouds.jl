#-
using Peridynamics #hide
using CairoMakie #src
CairoMakie.activate!(px_per_unit=2) #src
#-

# # [Point clouds](@id pointclouds)

# ## Point clouds
# In peridynamics, the continuum is mapped by material points.
# The point clouds are represented as [`PointCloud`](@ref) objects.
#
# ### Block with uniform distributed points:
# To generate a uniformly distributed `PointCloud` with the lengths `lx`, `ly`, `lz`, and
# point spacing `Δx`, simply type:
#-
lx1 = 3
ly1 = 1
lz1 = 1
Δx = 0.2
pc1 = PointCloud(lx1, ly1, lz1, Δx)
#-
fig = Figure(resolution = (1000,700), backgroundcolor = :transparent) #src
ax = Axis3(fig[1,1]; aspect =:data) #src
hidespines!(ax) #src
hidedecorations!(ax) #src
meshscatter!(ax, pc1.position; markersize=0.5Δx, color=:red) #src
save(joinpath(@__DIR__, "..", "assets", "pc1.png"), fig; px_per_unit=3) #src
fig #src
# ![](../assets/pc1.png) #md

#-
# The optional keyword arguments `center_x`, `center_y`, and `center_z` provide the
# possibility to specify the position of the [`PointCloud`](@ref) center.
# As seen in the image below, `pc2` (blue) is positioned with the `center`-keywords so that it forms a L-shape together with `pc1` (red).
#-
lx2 = 1
ly2 = 1
lz2 = 2
pc2 = PointCloud(lx2, ly2, lz2, Δx; center_x=(lx1-lx2)/2, center_z=(lz1+lz2)/2)
#-
fig = Figure(resolution = (700,700), backgroundcolor = :transparent) #src
ax = Axis3(fig[1,1]; aspect=:data) #src
hidespines!(ax) #src
hidedecorations!(ax) #src
meshscatter!(ax, pc1.position; markersize=0.5Δx, color=:red) #src
meshscatter!(ax, pc2.position; markersize=0.5Δx, color=:blue) #src
save(joinpath(@__DIR__, "..", "assets", "pc1_and_pc2.png"), fig; px_per_unit=3) #src
fig #src
# ![](../assets/pc1_and_pc2.png) #md
#-

# ### Merging of multiple point clouds

# The point clouds `pc1` and `pc2` can be merged to create one L-shaped `PointCloud` for a
# single body simulation. That can be accomplished with the [`pcmerge`](@ref) function:
#-
pc = pcmerge([pc1, pc2])
#- plot 3
fig = Figure(resolution = (700,700), backgroundcolor = :transparent) #src
ax = Axis3(fig[1,1]; aspect=:data) #src
hidespines!(ax) #src
hidedecorations!(ax) #src
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:green) #src
save(joinpath(@__DIR__, "..", "assets", "pc.png"), fig; px_per_unit=3) #src
fig #src
# ![](../assets/pc.png) #md
#-

# ### [Filtering of points regarding their position](@id filtering_points)
# To generate more complicated geometries from a uniform distributed block, points can be
# filtered out. For example, we want to model a cylinder with diameter $\text{\O}$ and
# thickness $t$:
#-
Ø = 1
t = 0.1
Δx = 0.03
pc0 = PointCloud(Ø, Ø, t, Δx)
#-
fig = Figure(backgroundcolor = :transparent, resolution = (1000, 500)) #src
ax = Axis3(fig[1,1]; aspect=:data) #src
hidespines!(ax) #src
hidedecorations!(ax) #src
meshscatter!(ax, pc0.position; markersize=0.5Δx, color=:blue) #src
save(joinpath(@__DIR__, "..", "assets", "pc0.png"), fig; px_per_unit=3) #src
fig #src
# ![](../assets/pc0.png) #md
#-

# Now we filter every point, that lies outside of the cylinder.
# Therefore, we search for all points that match the condition
# ```math
# \sqrt{{x_p}^2 + {y_p}^2} \leq \frac{\text{\O}}{2} \; ,
# ```
# with the $x$- and $y$-coordinate $x_p$ and $y_p$ of each point.
# The variable `cyl_id` contains the index of each point that matches this condition.
# Then we create a new point cloud `pc` using only the points of `pc0` specified in `cyl_id`.
#-
cyl_id = findall(p -> sqrt(p[1]^2 + p[2]^2) <= Ø/2, eachcol(pc0.position))
pc = PointCloud(pc0.position[:,cyl_id], pc0.volume[cyl_id])
#-
fig = Figure(backgroundcolor = :transparent, resolution = (1000, 500)) #src
ax = Axis3(fig[1,1]; aspect=:data) #src
hidespines!(ax) #src
hidedecorations!(ax) #src
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:blue) #src
save(joinpath(@__DIR__, "..", "assets", "pc0_filtered.png"), fig; px_per_unit=3) #src
fig #src
# ![](../assets/pc0_filtered.png) #md
#-
