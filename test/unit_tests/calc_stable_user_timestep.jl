using Peridynamics, Test

# TODO: test all materials and create the simmodel directly...
if Threads.nthreads() <= 2
    positions = [
        0.0 1.0
        0.0 0.0
        0.0 0.0
    ]
    Δx = 1.0
    δ = 1.5 * Δx
    n_points = 2
    volumes = fill(Δx^3, n_points)
    pc = PointCloud(positions, volumes)
    mat = BondBasedMaterial(; horizon=δ, rho=1, E=1, Gc=1)

    Δt = 0.7 * sqrt(2 * 1 / (Δx^3 * 1 / 1 * 18 * 2/3 / (π * δ^4)))
    @test calc_stable_user_timestep(pc, mat) == Δt
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
