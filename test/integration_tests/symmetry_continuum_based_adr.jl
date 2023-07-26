using Peridynamics, Test

# create perfect symmetrical PointCloud
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

# find the points that should have symmetrical positions after calculation
xp = findall(positions[1,:] .> 0)
xm = findall(positions[1,:] .< 0)
yp = findall(positions[2,:] .> 0)
ym = findall(positions[2,:] .< 0)
zp = findall(positions[3,:] .> 0)
zm = findall(positions[3,:] .< 0)

# define the material
mat = ContinuumBasedMaterial(horizon = 3.015point_spacing, rho=7850.0, E=210e9, nu=0.3, Gc=1.0)

# boundary conditions
fun1(t) = 0.0001
fun2(t) = -0.0001
bc1 = VelocityBC(fun1,findall(positions[3,:] .> width/2-0.6point_spacing),3)
bc2 = VelocityBC(fun2,findall(positions[3,:] .< -width/2+0.6point_spacing),3)
bcs = [bc1, bc2]

# settings
td = DynamicRelaxation(100)
es = ExportSettings(@__DIR__, 10)

# simulation job
job = PDSingleBodyAnalysis(name="MWE_Symmetry",pc=pc,mat=mat,bcs=bcs,td=td,es=es)

# remove all simulation files
rm.(filter(x->endswith(x,".vtu"), readdir(@__DIR__; join=true)), force=true)
rm.(filter(x->endswith(x,".log"), readdir(@__DIR__; join=true)), force=true)

# submit the simulation
model = submit(job)

# check if the correct files were exported
@test length(filter(x->endswith(x,".vtu"), readdir(@__DIR__))) == 11
@test length(filter(x->endswith(x,".log"), readdir(@__DIR__))) == 1

# check the symmetry
@test sort(model.position[1,xp]) ≈ sort(-model.position[1,xm])
@test sort(model.position[2,xp]) ≈ sort(model.position[2,xm])
@test sort(model.position[3,xp]) ≈ sort(model.position[3,xm])
@test sort(model.position[1,yp]) ≈ sort(model.position[1,ym])
@test sort(model.position[2,yp]) ≈ sort(-model.position[2,ym])
@test sort(model.position[3,yp]) ≈ sort(model.position[3,ym])
@test sort(model.position[1,zp]) ≈ sort(model.position[1,zm])
@test sort(model.position[2,zp]) ≈ sort(model.position[2,zm])
@test sort(model.position[3,zp]) ≈ sort(-model.position[3,zm])

# remove all simulation files
rm.(filter(x->endswith(x,".vtu"), readdir(@__DIR__; join=true)), force=true)
rm.(filter(x->endswith(x,".log"), readdir(@__DIR__; join=true)), force=true)
