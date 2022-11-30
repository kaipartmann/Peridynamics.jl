using Peridynamics
using Test

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

##------------------------------------------------------------------------------------------
# Contact Example with 4 points

if Threads.nthreads() <= 2
    positions1 = [
        0.0 1.0
        0.0 0.0
        0.0 0.0
    ]
    positions2 = [
        1.5 2.5
        0.0 0.0
        0.0 0.0
    ]
    point_spacing = 1.0
    δ = 1.5 * point_spacing
    E = 210e9
    rho = 7850.0
    n_points = 2
    volumes = fill(point_spacing^3, n_points)
    pc1 = PointCloud(positions1, volumes)
    pc2 = PointCloud(positions2, volumes)
    mat = BondBasedMaterial(horizon=δ, rho=rho, E=E, Gc=1.0)
    bodies = [BodySetup(pc1, mat), BodySetup(pc2, mat)]
    contact1 = Contact((1, 2), 1.0)
    job1 = PDContactAnalysis(;
        name="job1",
        body_setup=bodies,
        contact=[contact1],
        td=TimeDiscretization(1),
        es=ExportSettings(@__DIR__, 1),
    )
    rm.(joinpath.(@__DIR__,filter(x->endswith(x,".vtu"), readdir(@__DIR__))), force=true)
    rm.(joinpath.(@__DIR__,filter(x->endswith(x,".log"), readdir(@__DIR__))), force=true)
    bodies = submit(job1)

    @test contact1.spring_constant == 1e12
    b_int_contact = 0.5 * 9 * 1e12 / (π * δ^4) * (1-0.5) / 0.5
    @test bodies[1].b_int[1,2,1] == -b_int_contact
    @test bodies[2].b_int[1,1,1] == b_int_contact
    acc_contact = b_int_contact / rho
    @test bodies[1].acceleration[1,2,1] == -acc_contact
    @test bodies[2].acceleration[1,1,1] == acc_contact
    vel_contact = acc_contact * job1.td.Δt * 0.5
    @test bodies[1].velocity[1,2,1] == -vel_contact
    @test bodies[2].velocity[1,1,1] == vel_contact
    @test length(filter(x->endswith(x,".vtu"), readdir(@__DIR__))) == 4
    @test length(filter(x->endswith(x,".log"), readdir(@__DIR__))) == 1

    rm.(joinpath.(@__DIR__,filter(x->endswith(x,".vtu"), readdir(@__DIR__))), force=true)
    rm.(joinpath.(@__DIR__,filter(x->endswith(x,".log"), readdir(@__DIR__))), force=true)
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
