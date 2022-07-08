using Peridynamics
using Test

point_spacing = 0.25
width = 1
grid⁺ = point_spacing/2:point_spacing:width/2
grid⁻ = -reverse(grid⁺)
grid = [grid⁻; grid⁺]
positions = hcat(([x;y;z] for x in grid for y in grid for z in grid)...)
n_points = size(positions, 2)
volumes = fill(point_spacing^3, n_points)
pc = PointCloud(positions, volumes)
pc.failure_flag .= false # no failure possible in entire model!
xp = findall(positions[1,:] .> 0)
xm = findall(positions[1,:] .< 0)
yp = findall(positions[2,:] .> 0)
ym = findall(positions[2,:] .< 0)
zp = findall(positions[3,:] .> 0)
zm = findall(positions[3,:] .< 0)
mat = BondBasedMaterial(horizon = 3.015point_spacing, rho=7850.0, E=210e9, Gc=1.0)
fun1(t) = 100.
fun2(t) = -100.
bc1 = VelocityBC(fun1,findall(positions[3,:] .> width/2-0.6point_spacing),3)
bc2 = VelocityBC(fun2,findall(positions[3,:] .< -width/2+0.6point_spacing),3)
bcs = [bc1, bc2]
td = TimeDiscretization(100)
es = ExportSettings(@__DIR__, 10)
job = PDSingleBodyAnalysis(name="MWE_Symmetry",pc=pc,mat=mat,bcs=bcs,td=td,es=es)
rm.(joinpath.(@__DIR__,filter(x->endswith(x,".vtu"), readdir(@__DIR__))), force=true)
rm.(joinpath.(@__DIR__,filter(x->endswith(x,".jld2"), readdir(@__DIR__))), force=true)
rm.(joinpath.(@__DIR__,filter(x->endswith(x,".log"), readdir(@__DIR__))), force=true)
model = submit(job)

@test length(filter(x->endswith(x,".vtu"), readdir(@__DIR__))) == 11
@test length(filter(x->endswith(x,".jld2"), readdir(@__DIR__))) == 11
@test length(filter(x->endswith(x,".log"), readdir(@__DIR__))) == 1
@test sort(model.position[1,xp]) ≈ sort(-model.position[1,xm])
@test sort(model.position[2,xp]) ≈ sort(model.position[2,xm])
@test sort(model.position[3,xp]) ≈ sort(model.position[3,xm])
@test sort(model.position[1,yp]) ≈ sort(model.position[1,ym])
@test sort(model.position[2,yp]) ≈ sort(-model.position[2,ym])
@test sort(model.position[3,yp]) ≈ sort(model.position[3,ym])
@test sort(model.position[1,zp]) ≈ sort(model.position[1,zm])
@test sort(model.position[2,zp]) ≈ sort(model.position[2,zm])
@test sort(model.position[3,zp]) ≈ sort(-model.position[3,zm])

rm.(joinpath.(@__DIR__,filter(x->endswith(x,".vtu"), readdir(@__DIR__))), force=true)
rm.(joinpath.(@__DIR__,filter(x->endswith(x,".jld2"), readdir(@__DIR__))), force=true)
rm.(joinpath.(@__DIR__,filter(x->endswith(x,".log"), readdir(@__DIR__))), force=true)
