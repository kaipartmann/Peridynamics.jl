
#md # # [Spatial discretization](@id spatial_discretization)
#nb # # Spatial discretization

using Peridynamics #hide
using CairoMakie #hide
CairoMakie.activate!(px_per_unit=3) #hide

# ## Point clouds

# ### Block with uniform distributed points:

#- code block 1
lx1 = 3
ly1 = 1
lz1 = 1
Δx = 0.2
pc1 = PointCloud(lx1, ly1, lz1, Δx)

#- plot 1
fig = Figure(resolution = (1000,700), backgroundcolor = :transparent) #src
ax = Axis3(fig[1,1]; aspect =:data) #src
hidespines!(ax) #src
hidedecorations!(ax) #src
meshscatter!(ax, pc1.position; markersize=0.5Δx, color=:red) #src
save("docs/src/assets/pc1.png", fig; px_per_unit=3) #src
fig #src
#md # ![](../assets/pc1.png)

#- code block 2
lx2 = 1
ly2 = 1
lz2 = 2
pc2 = PointCloud(lx2, ly2, lz2, Δx; center_x=(lx1-lx2)/2, center_z=(lz1+lz2)/2)

#- plot 2
fig = Figure(backgroundcolor = :transparent)
ax = Axis3(fig[1,1]; aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
meshscatter!(ax, pc1.position; markersize=0.5Δx, color=:red)
meshscatter!(ax, pc2.position; markersize=0.5Δx, color=:blue)
fig

# ### Merging of multiple point clouds

#- code block 3
pc = pcmerge([pc1, pc2])

#- plot 3
fig = Figure(backgroundcolor = :transparent)
ax = Axis3(fig[1,1]; aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:green)
fig

# ### Filtering of points regarding their position

#- code block 4
Ø = 1
t = 0.1
Δx = 0.03
pc0 = PointCloud(Ø, Ø, t, Δx)

#-
fig = Figure(backgroundcolor = :transparent, resolution = (1000, 500))
ax = Axis3(fig[1,1]; aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
meshscatter!(ax, pc0.position; markersize=0.5Δx, color=:blue)
fig

#-
cyl_id = findall(p -> sqrt(p[1]^2 + p[2]^2) <= Ø/2, eachcol(pc0.position))
pc = PointCloud(pc0.position[:,cyl_id], pc0.volume[cyl_id])

#-
fig = Figure(backgroundcolor = :transparent, resolution = (1000, 500))
ax = Axis3(fig[1,1]; aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:blue)
fig
