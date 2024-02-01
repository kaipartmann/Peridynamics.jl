## single body simulation

# p::DiscretizationSetup
#    -> PointCloud
#        -> generate with lx, ly, lz, ΔX0
#        -> generate with point position & volume
#    -> Mesh (Tetraeder & Hexaeder)
#        -> read from Abaqus-file -> generate Points and volume from mesh
#        ->
p = PointCloud(...)

body = Body(p::DiscretizationSetup)
















# old
using Peridynamics

lx, lyz, Δx = 0.2, 0.002, 0.002/4
pc = PointCloud(lx, lyz, lyz, Δx)
pc.failure_flag .= false

mat = BondBasedMaterial(horizon=3.015Δx, rho=7850, E=210e9, epsilon_c=0.01)

T, vmax = 1.0e-5, 2
points_fix_left = findall(pc.position[1,:] .< -lx/2+1.2Δx)
vwave(t) = t < T ? vmax * sin(2π/T * t) : 0
bcs = [VelocityBC(vwave, points_fix_left, 1)]

es = ExportSettings(resfolder, 10)
vv = VelocityVerlet(2000)
job = PDSingleBodyAnalysis(name=simname, pc=pc, mat=mat, bcs=bcs, es=es, td=vv)
submit(job)

















# NEW API
using Peridynamics

lx, lyz, Δx = 0.2, 0.002, 0.002/4
body = Body(lx, lyz, lyz, Δx)
failure_allowed!(body, false)
material!(body, BBMaterial(horizon=3.015Δx, rho=7850, E=210e9, epsilon_c=0.01))
point_set!((x,y,z) -> x < -lx/2+1.2Δx, body, :fix_left)
failure_allowed!(body, :fix_left, false)
velocity_bc!(t -> t < 1e-5 ? 2 * sin(2π/1e-5 * t) : 0, body, :fix_left, :x)
solver = VelocityVerlet(steps=2000)
config = ExportSettings(path="results", freq=10)
job = Job(body, solver, config)
submit(job)

# NEW API CONTACT
using Peridynamics
lx, lyz, Δx = 0.2, 0.002, 0.002/4
b1 = Body(lx, lyz, lyz, Δx)
b2 = Body(lx, 2lyz, 2lyz, 2Δx)
failure_allowed!(b1, false)
material!(b1, BBMaterial(horizon=3.015Δx, rho=7850, E=210e9, epsilon_c=0.01))
material!(b2, BBMaterial(horizon=3.015Δx, rho=7850, E=210e9, epsilon_c=0.01))
point_set!(p -> p[1] < -lx/2+1.2Δx, b1, :fix_left)
velocity_ic!(body, :fix_left, 1, -100)
vv = VelocityVerlet(steps=2000)
es = ExportSettings(path="results", freq=10)
job = Job(b1, b2; solver=VelocityVerlet(steps=2000), config=Settings(path="results", freq=10))
submit(job)







####
using Peridynamics

b = Body(rand(3,10), rand(10))

point_set!(b, :a, 1:3)

b.psh.point_sets

mat = BBMaterial(horizon=1, E=1, rho=1, Gc=1)

material!(b, mat)

velocity_bc!(t -> t, b, :a, 2)
forcedensity_bc!(t -> t, b, :a, 1)

b.conditions

##
function testargs(a)

    return Dict(a)
end



##

point_set!(x -> x > 0, b, :set_a)
    points = findall(x -> x > 0, b.position[1,:])
