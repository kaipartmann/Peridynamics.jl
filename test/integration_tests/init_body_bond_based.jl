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
    body = Peridynamics.init_body(mat, pc)

    # show BondBasedBody
    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (20, 40)), "text/plain", body)
    msg_body = String(take!(io))
    @test msg_body == "2-points Peridynamics.BondBasedBody:\n  1 bonds"

    #TODO: test body!!!
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
