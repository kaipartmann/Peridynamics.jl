using Peridynamics, Test

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
    mat = BBMaterial(horizon=δ, rho=rho, E=E, Gc=1.0)
    bodies = [BodySetup(pc1, mat), BodySetup(pc2, mat)]
    contact1 = Contact((1, 2), 1.0)
    job1 = PDContactAnalysis(;
        name="job1",
        body_setup=bodies,
        contact=[contact1],
        td=VelocityVerlet(1),
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
