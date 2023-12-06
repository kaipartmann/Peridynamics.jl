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
    mat = CPDMaterial(horizon=δ, rho=7850.0, E=E, nu=0.3, Gc=1.0)
    body = Peridynamics.init_body(mat, pc)

    # show BondBasedBody
    io = IOBuffer()
    show(io, "text/plain", body)
    msg_body = String(take!(io))
    @test msg_body == "2-points Peridynamics.ContinuumBasedSimBody:\n  2 bonds\n  0 two-neighbor interactions\n  0 three-neighbor interactions\n"

    #TODO: test body!!!
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
