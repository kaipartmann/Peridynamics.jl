using Peridynamics
using Test

##------------------------------------------------------------------------------------------
# PointCloud

position1 = [
    1.0 2.0
    3.0 4.0
    5.0 6.0
]
volume1 = [1.0, 1.0]
pc1 = PointCloud(position1, volume1)
@test pc1.n_points == 2
@test pc1.failure_flag == [true, true]
@test pc1.radius == Peridynamics.sphere_radius.(volume1)
@test isempty(pc1.point_sets)
@test_throws DimensionMismatch PointCloud([1.0 2.0; 3.0 4.0], volume1)

pc2 = PointCloud(1, 1, 1, 0.5)
@test pc2.position == [
    -0.25  -0.25  -0.25  -0.25   0.25   0.25   0.25  0.25
    -0.25  -0.25   0.25   0.25  -0.25  -0.25   0.25  0.25
    -0.25   0.25  -0.25   0.25  -0.25   0.25  -0.25  0.25
]
@test pc2.n_points == 8
@test pc2.volume == fill(0.125, 8)
@test pc2.radius == Peridynamics.sphere_radius.(fill(0.125, 8))
@test pc2.failure_flag == fill(true, 8)
@test isempty(pc2.point_sets)

pc3 = PointCloud(1, 1, 1, 0.5; center_x = 1, center_y = 1, center_z = 1)
@test pc3.position == [
    0.75  0.75  0.75  0.75  1.25  1.25  1.25  1.25
    0.75  0.75  1.25  1.25  0.75  0.75  1.25  1.25
    0.75  1.25  0.75  1.25  0.75  1.25  0.75  1.25
]
@test pc3.n_points == 8
@test pc3.volume == fill(0.125, 8)
@test pc3.radius == Peridynamics.sphere_radius.(fill(0.125, 8))
@test pc3.failure_flag == fill(true, 8)
@test isempty(pc3.point_sets)

##------------------------------------------------------------------------------------------
# TimeDiscretization

td1 = TimeDiscretization(1)
@test td1.n_timesteps == 1
@test td1.Δt == -1
@test td1.alg == :verlet

td2 = TimeDiscretization(2; alg=:dynrelax)
@test td2.n_timesteps == 2
@test td2.Δt == 1
@test td2.alg == :dynrelax

@test_throws DomainError TimeDiscretization(2; alg = :somethingwrong)

td3 = TimeDiscretization(3, 0.1)
@test td3.n_timesteps == 3
@test td3.Δt == 0.1
@test td3.alg == :verlet

msg = "for dynamic relaxation a time step of Δt = 1 is recommended!"
td4 = @test_warn msg TimeDiscretization(3, 0.2; alg=:dynrelax)
@test td4.n_timesteps == 3
@test td4.Δt == 0.2
@test td4.alg == :dynrelax

@test_throws DomainError TimeDiscretization(4, 0.4; alg = :somethingwrong)


##------------------------------------------------------------------------------------------
# ExportSettings

es1 = ExportSettings()
@test es1.path == ""
@test es1.exportfreq == 0
@test es1.resultfile_prefix == ""
@test es1.logfile == ""
@test es1.exportflag == false

es2 = ExportSettings("test/path", 10)
@test es2.path == "test/path"
@test es2.exportfreq == 10
@test es2.resultfile_prefix == ""
@test es2.logfile == ""
@test es2.exportflag == true

##------------------------------------------------------------------------------------------
# Boundary Conditions and Initial Conditions

vbc1 = VelocityBC(t -> 1.0, 1:5, 1)
@test vbc1.fun(0) == 1.0
@test_throws MethodError vbc1.fun() == 1.0
@test vbc1.point_id_set == [1, 2, 3, 4, 5]
@test vbc1.dim == 1

vic1 = VelocityIC(1.0, 1:5, 1)
@test vic1.val == 1.0
@test vic1.point_id_set == [1, 2, 3, 4, 5]
@test vic1.dim == 1

fdbc1 = ForceDensityBC(t -> 1.0, 1:5, 1)
@test fdbc1.fun(0) == 1.0
@test_throws MethodError fdbc1.fun() == 1.0
@test fdbc1.point_id_set == [1, 2, 3, 4, 5]
@test fdbc1.dim == 1

pdvbc1 = PosDepVelBC((x, y, z, t) -> 1.0, 1:5, 1)
@test pdvbc1.fun(0, 0, 0, 0) == 1.0
@test_throws MethodError pdvbc1.fun() == 1.0
@test pdvbc1.point_id_set == [1, 2, 3, 4, 5]
@test pdvbc1.dim == 1

##------------------------------------------------------------------------------------------
# Contact

c1 = Contact((1,2), 0.1)
@test c1.body_id_set == (1,2)
@test c1.search_radius == 0.1
@test c1.spring_constant == 1e12

c2 = Contact((3,4), 0.2, 2e13)
@test c2.body_id_set == (3,4)
@test c2.search_radius == 0.2
@test c2.spring_constant == 2e13

##------------------------------------------------------------------------------------------
# BodySetup

bs1 = BodySetup(PointCloud(1, 1, 1, 0.5), BondBasedMaterial(; horizon=1, rho=1, E=1, Gc=1))
@test isempty(bs1.precracks)
@test isempty(bs1.bcs)
@test isempty(bs1.ics)

##------------------------------------------------------------------------------------------
# Peridynamic Analysis structs

pdsba1 = PDSingleBodyAnalysis(;
    name="pdsba1",
    pc=PointCloud(1, 1, 1, 0.5),
    mat=BondBasedMaterial(; horizon=1, rho=1, E=1, Gc=1),
    td=TimeDiscretization(10),
    es=ExportSettings(),
)
@test pdsba1.es.resultfile_prefix == "pdsba1"
@test pdsba1.es.logfile =="pdsba1.log"
@test pdsba1.es.exportfreq == 11
@test isempty(pdsba1.precracks)
@test isempty(pdsba1.bcs)
@test isempty(pdsba1.ics)

pdsba2 = PDSingleBodyAnalysis(;
    name="pdsba2",
    pc=PointCloud(1, 1, 1, 0.5),
    mat=BondBasedMaterial(; horizon=1, rho=1, E=1, Gc=1),
    td=TimeDiscretization(10),
    es=ExportSettings("test",2),
)
@test pdsba2.es.resultfile_prefix == joinpath("test", "pdsba2")
@test pdsba2.es.logfile == joinpath("test", "pdsba2.log")
@test pdsba2.es.exportfreq == 2
@test isempty(pdsba2.precracks)
@test isempty(pdsba2.bcs)
@test isempty(pdsba2.ics)

pdca1 = PDContactAnalysis(;
    name="pdca1",
    body_setup=[bs1, bs1],
    contact=[c1],
    td=TimeDiscretization(10),
    es=ExportSettings(),
)
@test pdca1.es.resultfile_prefix == "pdca1"
@test pdca1.es.logfile == "pdca1.log"
@test pdca1.es.exportfreq == 11

pdca2 = PDContactAnalysis(;
    name="pdca2",
    body_setup=[bs1, bs1],
    contact=[c1],
    td=TimeDiscretization(10),
    es=ExportSettings("test",2),
)
@test pdca2.es.resultfile_prefix == joinpath("test", "pdca2")
@test pdca2.es.logfile == joinpath("test", "pdca2.log")
@test pdca2.es.exportfreq == 2

@test_throws ErrorException PDContactAnalysis(;
    name="pdca2",
    body_setup=[bs1, bs1],
    contact=[c1],
    td=TimeDiscretization(10; alg=:dynrelax),
    es=ExportSettings("test",2),
)
