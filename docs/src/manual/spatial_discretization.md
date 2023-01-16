# [Spatial discretization](@id spatial_discretization)
### Point clouds

In peridynamics, the continuum is mapped by material points.
The point clouds are represented as [`PointCloud`](@ref) objects.

##### Uniform distributed block of points
To generate a uniformly distributed `PointCloud` with the lengths `lx`, `ly`, `lz`, and
point spacing `Δx`, simply type:
```@setup spatial_discretization_1
import Pkg
Pkg.activate(joinpath("..", "..","Project.toml"))
using Peridynamics
```
```@example spatial_discretization_1
lx = 3
ly = 1
lz = 1
Δx = 0.2
pc = PointCloud(lx, ly, lz, Δx)
```

The points can be easily displayed with [Makie.jl](https://docs.makie.org/stable/):
```@setup spatial_discretization_1
using WGLMakie
using JSServe
Page(exportable=true, offline=true)
WGLMakie.activate!()
```
```@example spatial_discretization_1
fig = Figure() #hide
ax = Axis3(fig[1,1]; aspect=:data) #hide
hidespines!(ax) #hide
hidedecorations!(ax) #hide
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:red) #hide
fig #hide
```

The optional keyword arguments `center_x`, `center_y`, and `center_z` provide the possibility to specify the position of the [`PointCloud`](@ref) center.
```@example spatial_discretization_1
pc2 = PointCloud(1, 1, 2, Δx; center_x=1, center_z=1.5)
```
```@example spatial_discretization_1
fig = Figure() #hide
ax = Axis3(fig[1,1]; aspect=:data) #hide
hidespines!(ax) #hide
hidedecorations!(ax) #hide
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:red) #hide
meshscatter!(ax, pc2.position; markersize=0.5Δx, color=:blue) #hide
fig #hide
```


### Cracks & damage
### Read Abaqus mesh