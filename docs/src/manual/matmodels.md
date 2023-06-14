# [Material models](@id mat_models)

## MultiMaterial

```@setup mm1
using Peridynamics
using CairoMakie
```
```@example mm1
lx = 3
ly = 1
lz = 1
Δx = 0.15
pc = PointCloud(lx, ly, lz, Δx)
pc.failure_flag .= false
```

```@example mm1
fig = Figure()
ax = Axis3(fig[1,1]; aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
meshscatter!(ax, pc.position; markersize=0.5Δx, color=:blue)
fig
```

```@example mm1
mat1 = BondBasedMaterial(horizon=3.1Δx, rho=2000, E=20e9, Gc=10)
```

```@example mm1
mat2 = BondBasedMaterial(horizon=3.1Δx, rho=7800, E=210e9, Gc=600)
```

```@example mm1
mat_of_point = ifelse.(pc.position[1,:] .> 0, 1, 2)
```


```@example mm1
mm = MultiMaterial((mat1, mat2), mat_of_point)
```

```@example mm1
id_mat1 = findall(x -> x == 1, mat_of_point)
id_mat2 = findall(x -> x == 2, mat_of_point)
```

```@example mm1
fig = Figure()
ax = Axis3(fig[1,1]; aspect=:data)
hidespines!(ax)
hidedecorations!(ax)
pltlg1 = scatter!(ax, pc.position[:, id_mat1[1]]; color=:blue, label="mat1")
pltlg2 = scatter!(ax, pc.position[:, id_mat2[1]]; color=:red, label="mat2")
meshscatter!(ax, pc.position[:, id_mat1]; markersize=0.5Δx, color=:blue)
meshscatter!(ax, pc.position[:, id_mat2]; markersize=0.5Δx, color=:red)

axislegend(ax; position=:lt)

pltlg1.visible[] = false
pltlg2.visible[] = false
fig
```

```@example mm1
id_left = findall(pc.position[1,:] .> lx/2-Δx)
id_right = findall(pc.position[1,:] .< -lx/2+Δx)

bcs = [
    ForceDensityBC(t -> 1e8, id_left, 1),
    ForceDensityBC(t -> 1e8, id_right, 1),
]

job = PDSingleBodyAnalysis(;
    name="test",
    pc=pc,
    mat=mm,
    bcs=bcs,
    td=TimeDiscretization(500),
    es=ExportSettings(),
)

resbody = submit(job);
```