# [Spatial discretization](@id spatial_discretization)
## Point clouds

In peridynamics, the continuum is mapped by material points.
The point clouds are represented as [`PointCloud`](@ref) objects.

### Uniform distributed block of points
To generate a uniformly distributed `PointCloud` with the lengths `lx`, `ly`, `lz`, and
point spacing `Δx`, simply type:
```@setup sd1
using Peridynamics
using CairoMakie
```
```@example sd1
lx1 = 3
ly1 = 1
lz1 = 1
Δx = 0.2
pc1 = PointCloud(lx1, ly1, lz1, Δx)
```
The points can be easily displayed with [Makie.jl](https://docs.makie.org/stable/):
```@raw html
<details>
    <summary>Code for the plot:</summary>
```
```julia
fig = Figure()
ax = Axis3(fig[1,1]; aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
meshscatter!(ax, pc1.position; markersize=0.5Δx, color=:red)
fig
```
```@raw html
</details>
```
```@example sd1
fig = Figure() #hide
ax = Axis3(fig[1,1]; aspect=:data) #hide
hidespines!(ax) #hide
hidedecorations!(ax) #hide
meshscatter!(ax, pc1.position; markersize=0.5Δx, color=:red) #hide
fig #hide
```

The optional keyword arguments `center_x`, `center_y`, and `center_z` provide the possibility to specify the position of the [`PointCloud`](@ref) center.
As seen in the image below, `pc2` is positioned with the `center`-keywords so that it forms a L-shape together with `pc1`.
```@example sd1
lx2 = 1
ly2 = 1
lz2 = 2
pc2 = PointCloud(lx2, ly2, lz2, Δx; center_x=(lx1-lx2)/2, center_z=(lz1+lz2)/2)
```
```@raw html
<details>
    <summary>Code for the plot:</summary>
```
```julia
fig = Figure()
ax = Axis3(fig[1,1]; aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
meshscatter!(ax, pc1.position; markersize=0.5Δx, color=:red)
meshscatter!(ax, pc2.position; markersize=0.5Δx, color=:blue)
fig
```
```@raw html
</details>
```
```@example sd1
fig = Figure() #hide
ax = Axis3(fig[1,1]; aspect=:data) #hide
hidespines!(ax) #hide
hidedecorations!(ax) #hide
meshscatter!(ax, pc1.position; markersize=0.5Δx, color=:red) #hide
meshscatter!(ax, pc2.position; markersize=0.5Δx, color=:blue) #hide
fig #hide
```

### Merging of multiple point clouds
The point clouds `pc1` and `pc2` can be merged to create one L-shaped `PointCloud` for a single body simulation.
That can be accomplished with the [`pcmerge`](@ref) function:
```@example sd1
pc = pcmerge([pc1, pc2])
```
```@raw html
<details>
    <summary>Code for the plot:</summary>
```
```julia
fig = Figure()
ax = Axis3(fig[1,1]; aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:green)
fig
```
```@raw html
</details>
```
```@example sd1
fig = Figure() #hide
ax = Axis3(fig[1,1]; aspect=:data) #hide
hidespines!(ax) #hide
hidedecorations!(ax) #hide
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:green) #hide
fig #hide
```

### Filter points regarding their position
To generate more complicated geometries from a uniform distributed block, points can be filtered out.
For example, we want to model a cylinder with diameter $\text{\O}$ and thickness $t$.
```@example sd1
Ø = 1
t = 0.1
Δx = 0.03
pc0 = PointCloud(Ø, Ø, t, Δx)
```
```@raw html
<details>
    <summary>Code for the plot:</summary>
```
```julia
fig = Figure()
ax = Axis3(fig[1,1]; aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
meshscatter!(ax, pc0.position; markersize=0.5Δx, color=:blue)
display(fig)
```
```@raw html
</details>
```
```@example sd1
fig = Figure() #hide
ax = Axis3(fig[1,1]; aspect=:data) #hide
hidespines!(ax) #hide
hidedecorations!(ax) #hide
meshscatter!(ax, pc0.position; markersize=0.5Δx, color=:blue) #hide
fig #hide
```

Now we filter every point, that lies outside of the cylinder.
Therefore, we search for all points that match the condition
```math
\sqrt{{x_p}^2 + {y_p}^2} \leq \frac{\text{\O}}{2} \; ,
```
with the $x$- and $y$-coordinate $x_p$ and $y_p$ of each point.
The variable `cyl_id` contains the index of each point that matches this condition.
Then we create a new point cloud `pc` using only the points of `pc0` specified in `cyl_id`.
```@example sd1
cyl_id = findall(p -> sqrt(p[1]^2 + p[2]^2) <= Ø/2, eachcol(pc0.position))
pc = PointCloud(pc0.position[:,cyl_id], pc0.volume[cyl_id])
```
```@raw html
<details>
    <summary>Code for the plot:</summary>
```
```julia
fig = Figure()
ax = Axis3(fig[1,1]; aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:blue)
fig
```
```@raw html
</details>
```
```@example sd1
fig = Figure() #hide
ax = Axis3(fig[1,1]; aspect=:data) #hide
hidespines!(ax) #hide
hidedecorations!(ax) #hide
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:blue) #hide
fig #hide
```
    
## Cracks & damage



## Read Abaqus mesh