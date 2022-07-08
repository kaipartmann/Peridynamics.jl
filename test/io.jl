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
    rm.(joinpath.(@__DIR__,filter(x->endswith(x,".vtu"), readdir(@__DIR__))), force=true)
    rm.(joinpath.(@__DIR__,filter(x->endswith(x,".jld2"), readdir(@__DIR__))), force=true)
    Peridynamics.export_results(body, "testfile", 0, 0.0)

    @test isfile("testfile_t0.vtu")
    @test isfile("testfile_t0.jld2")

    rm.(joinpath.(@__DIR__,filter(x->endswith(x,".vtu"), readdir(@__DIR__))), force=true)
    rm.(joinpath.(@__DIR__,filter(x->endswith(x,".jld2"), readdir(@__DIR__))), force=true)

    msg_mat = """
    BondBasedMaterial:
      δ: 1.5
      rho: 7850
      E: 2.1e+11
      nu: 0.25
      G: 8.4e+10
      K: 1.4e+11
      bc: 1.58448e+11
      Gc: 1
      εc: 1.6265e-06
    """
    msg_body = """
    2-points Peridynamics.BondBasedBody:
      Number of bonds / 1-NI: 1
    """
    @test_warn msg_mat Base.show(stderr, "text/plain", mat)
    @test_warn msg_body Base.show(stderr, "text/plain", body)
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
