@testitem "get_storage" begin
    position, volume = uniform_box(1,1,1,0.25)
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    ts = VelocityVerlet(steps=10)
    pd = Peridynamics.PointDecomposition(body, 1)
    system = Peridynamics.get_system(body, pd, 1)
    storage = Peridynamics.get_storage(body.mat, ts, system)

    @test storage isa Peridynamics.BBStorage
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
end

@testitem "material! BondBased" begin
    position, volume = uniform_box(1,1,1,0.5)
    body = Body(BBMaterial(), position, volume)
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

    material!(body, horizon=2, rho=1, lambda=0.4, Gc=1)
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

    material!(body, horizon=2, rho=1, mu=0.4, Gc=1)
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

    material!(body, horizon=2, rho=1, E=1, G=0.4, Gc=1)
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

    material!(body, horizon=2, rho=1, G=0.4, nu=0.25, Gc=1)
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

    material!(body, horizon=2, rho=1, lambda=0.4, mu=0.4, Gc=1)
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
end
