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

    @test body.n_points == n_points
    @test body.n_bonds == 1
    @test body.n_threads == Threads.nthreads()
    @test body.unique_bonds == true
    @test body.b_int == zeros(3, n_points, Threads.nthreads())
    if Threads.nthreads() == 1
        @test body.n_active_family_members == [1; 1;;]
    elseif Threads.nthreads() == 2
        @test body.n_active_family_members == [1 0; 0 1]
    end
    @test body.bond_failure == [1]

    # Boundary Condition:
    # Point 2 with ẋ_z = 1 m/s with Δt = 0.0015 s
    body.position[1, 2] = 1.0015
    Peridynamics.compute_forcedensity!(body, mat)

    if Threads.nthreads() == 1
        @test body.n_active_family_members == [0; 0;;]
    elseif Threads.nthreads() == 2
        @test body.n_active_family_members == [0 0; 0 0]
    end
    @test body.bond_failure == [0]

    b¹² = [
        18 * E / (3 * (1 - 2 * 0.25)) / (π * δ^4) * 1.0015 * 0.0015/1.0015 * 1.0
        0.0
        0.0
    ]
    Peridynamics.update_thread_cache!(body)
    @test body.b_int[:,:,1] ≈ [b¹² -b¹²]
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
