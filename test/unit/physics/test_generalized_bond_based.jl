# `GBBMaterial`: the generalized bond-based formulation of `src/physics/generalized_bond_based.jl`.

@testitem "get_storage GBBMaterial" begin
    position, volume = uniform_box(1,1,1,0.25)
    body = Body(GBBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    ts = VelocityVerlet(steps=10)
    pd = Peridynamics.PointDecomposition(body, 1)
    system = Peridynamics.get_system(body, pd, 1)
    storage = Peridynamics.get_storage(body.mat, ts, system)

    @test storage isa Peridynamics.GBBStorage
    @test storage.position == position
    @test storage.displacement == zeros(3, 64)
    @test storage.velocity == zeros(3, 64)
    @test storage.velocity_half == zeros(3, 64)
    @test storage.acceleration == zeros(3, 64)
    @test storage.b_int == zeros(3, 64)
    @test storage.b_ext == zeros(3, 64)
    @test storage.damage == zeros(64)
    @test storage.bond_active == ones(Bool, 4032)
    @test storage.n_active_bonds == fill(63, 64)
    @test storage.weighted_volume == zeros(64)
end

@testitem "material! GBBMaterial" begin
    position, volume = uniform_box(1,1,1,0.5)
    body = Body(GBBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1, nu=0.25, Gc=1)

    (; δ, rho, E, nu, G, K, λ, μ, Gc, εc) = body.point_params[1]
    @test δ ≈ 2.0
    @test rho ≈ 1.0
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 0.6666666666666666
    @test λ ≈ 0.4
    @test μ ≈ 0.4
    @test Gc ≈ 1.0
    @test εc ≈ 0.6454972243679028

    @test_throws ArgumentError material!(body, horizon=2, rho=1, E=1, nu=0.26, Gc=1)

    material!(body, horizon=2, rho=1, E=1, Gc=1)
    (; δ, rho, E, nu, G, K, λ, μ, Gc, εc) = body.point_params[1]
    @test δ ≈ 2.0
    @test rho ≈ 1.0
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 0.6666666666666666
    @test λ ≈ 0.4
    @test μ ≈ 0.4
    @test Gc ≈ 1.0
    @test εc ≈ 0.6454972243679028

    @test_throws ArgumentError material!(body, horizon=2, rho=1, nu=0.25, Gc=1)

    material!(body, horizon=2, rho=1, G=0.4, Gc=1)
    (; δ, rho, E, nu, G, K, λ, μ, Gc, εc) = body.point_params[1]
    @test δ ≈ 2.0
    @test rho ≈ 1.0
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 0.6666666666666666
    @test λ ≈ 0.4
    @test μ ≈ 0.4
    @test Gc ≈ 1.0
    @test εc ≈ 0.6454972243679028

    material!(body, horizon=2, rho=1, K=2/3, Gc=1)
    (; δ, rho, E, nu, G, K, λ, μ, Gc, εc) = body.point_params[1]
    @test δ ≈ 2.0
    @test rho ≈ 1.0
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 0.6666666666666666
    @test λ ≈ 0.4
    @test μ ≈ 0.4
    @test Gc ≈ 1.0
    @test εc ≈ 0.6454972243679028

    # parameters that lead to nu ≠ 0.25
    @test_throws ArgumentError material!(body, horizon=2, rho=1, E=1, G=1, Gc=1)
end

@testitem "GBBMaterial: force density of a stretched bond" begin
    ref_position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    δ = 1.5
    E = 210e9
    body = Body(GBBMaterial(), ref_position, volume)
    material!(body; horizon=δ, rho=1, E)
    no_failure!(body)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; position, b_int) = chunk.storage
    (; bonds) = chunk.system

    @test position == ref_position
    @test b_int == zeros(3, 2)

    # Boundary Condition:
    # Point 2 with v_z = 1 m/s with Δt = 0.0015 s
    position[1, 2] = 1.0015

    Peridynamics.calc_force_density!(chunk, 0, 0)

    wvol = 1.0
    b12 = 18 * E / (3 * (1 - 2 * 0.25)) / wvol * 1.0015 * 0.0015/1.0015 * 1.0
    @test b_int ≈ [b12 -b12; 0.0 0.0; 0.0 0.0]
end

@testitem "GBBMaterial: force density with energy surface correction" begin
    ref_position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    δ = 1.5
    E = 210e9
    body = Body(GBBMaterial{EnergySurfaceCorrection}(), ref_position, volume)
    material!(body; horizon=δ, rho=1, E)
    no_failure!(body)

    ts = VelocityVerlet(steps=1)
    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; position, b_int) = chunk.storage
    (; bonds) = chunk.system

    @test position == ref_position
    @test b_int == zeros(3, 2)

    # Boundary Condition:
    # Point 2 with v_z = 1 m/s with Δt = 0.0015 s
    position[1, 2] = 1.0015

    Peridynamics.calc_force_density!(chunk, 0, 0)

    sc = 0.2 # surface correction factor
    @test all(x -> x ≈ sc, chunk.system.correction.scfactor)

    wvol = 0.2
    b12 = sc * 18 * E / (3 * (1 - 2 * 0.25)) / wvol * 1.0015 * 0.0015/1.0015 * 1.0
    @test b_int ≈ [b12 -b12; 0.0 0.0; 0.0 0.0]
end

@testitem "GBBMaterial: force density across a material interface" begin
    ref_position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    δ = 1.5
    Ea = 105e9
    Eb = 2 * Ea
    E = (Ea + Eb) / 2
    body = Body(GBBMaterial(), ref_position, volume)
    point_set!(body, :a, [1])
    point_set!(body, :b, [2])
    material!(body, :a, horizon=δ, rho=1, E=Ea, Gc=1.0)
    material!(body, :b, horizon=δ, rho=1, E=Eb, Gc=1.0)
    no_failure!(body)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; position, b_int) = chunk.storage
    (; bonds) = chunk.system

    @test position == ref_position
    @test b_int == zeros(3, 2)

    # Boundary Condition:
    # Point 2 with v_z = 1 m/s with Δt = 0.0015 s
    position[1, 2] = 1.0015

    Peridynamics.calc_force_density!(chunk, 0, 0)

    wvol = 1.0
    b12 = 18 * E / (3 * (1 - 2 * 0.25)) / wvol * 1.0015 * 0.0015/1.0015 * 1.0
    @test b_int ≈ [b12 -b12; 0.0 0.0; 0.0 0.0]
end
