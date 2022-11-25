using Peridynamics
using Test

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
    Peridynamics.export_vtk(body, "testfile", 0, 0.0)

    @test isfile("testfile_t0.vtu")

    rm.(joinpath.(@__DIR__,filter(x->endswith(x,".vtu"), readdir(@__DIR__))), force=true)

    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (20, 40)), "text/plain", mat)
    msg_mat = String(take!(io))
    @test msg_mat == string(
        "BondBasedMaterial:\n  δ:   1.5\n  rho: 7850.0\n  E:   2.1e11\n  ",
        "nu:  0.25\n  G:   8.4e10\n  K:   1.4e11\n  bc:  1.584475877892647e11\n  ",
        "Gc:  1.0\n  εc:  1.6265001215808886e-6",
    )

    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (20, 40)), "text/plain", body)
    msg_body = String(take!(io))
    @test msg_body == "2-points Peridynamics.BondBasedBody:\n  1 bonds"
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
