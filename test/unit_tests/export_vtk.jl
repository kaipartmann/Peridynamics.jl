using Peridynamics, Test

if Threads.nthreads() <= 2
    positions = [
        0.0 1.0
        0.0 0.0
        0.0 0.0
    ]
    point_spacing = 1.0
    δ = 1.5 * point_spacing
    E = 210e9
    n_points = 2
    volumes = fill(point_spacing^3, n_points)
    pc = PointCloud(positions, volumes)
    mat = BondBasedMaterial(horizon=δ, rho=7850.0, E=E, Gc=1.0)
    body = Peridynamics.create_simmodel(mat, pc)

    # export vtk
    rm.(filter(x->endswith(x,".vtu"), readdir(@__DIR__; join=true)), force=true)
    Peridynamics.export_vtk(body, joinpath(@__DIR__, "testfile"), 0, 0.0)

    filename = joinpath(@__DIR__,"testfile_t0000.vtu")
    @test isfile(filename)

    results = read_vtk(filename)
    @test results["Position"] == positions
    @test results["Damage"] == zeros(2)
    @test results["Displacement"] == zeros(3,2)
    @test results["Velocity"] == zeros(3,2)
    @test results["ForceDensity"] == zeros(3,2)
    @test results["Time"] == zeros(1)

    rm.(filter(x->endswith(x,".vtu"), readdir(@__DIR__; join=true)), force=true)
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
