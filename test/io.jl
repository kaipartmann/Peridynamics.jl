using Peridynamics
using Test
using FileIO

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
    Peridynamics.export_results(body, "testfile", 0, 0.0)

    @test isfile("testfile_t0.vtu")
    @test isfile("testfile_t0.jld2")
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
