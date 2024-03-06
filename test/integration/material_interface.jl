
@testitem "TestMaterial halo exchange" begin
    include("test_material.jl")

    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = TestMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, nu=0.25, Gc=1)
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    velocity_ic!(body, :a, :x, 1.0)
    velocity_bc!(t -> t, body, :a, :x)
    forcedensity_bc!(t -> t, body, :a, :x)
    precrack!(body, :a, :b)
    ts = VelocityVerlet(steps=10)
    point_decomp = Peridynamics.PointDecomposition(body, 2)
    tdh = Peridynamics.ThreadsDataHandler(body, ts, point_decomp)

    b1 = tdh.chunks[1]
    @test b1 isa Peridynamics.BodyChunk
    @test b1.store.position == position
    @test b1.store.displacement == zeros(3, 2)
    @test b1.store.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.store.velocity_half == zeros(3, 2)
    @test b1.store.acceleration == zeros(3, 2)
    @test b1.store.b_int == zeros(3, 4)
    @test b1.store.b_ext == zeros(3, 2)
    @test b1.store.damage ≈ [2/3, 2/3]
    @test b1.store.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.store.n_active_bonds == [1, 1]
    b2 = tdh.chunks[2]
    @test b2 isa Peridynamics.BodyChunk
    @test b2.store.position == position[:, [3, 4, 1, 2]]
    @test b2.store.displacement == zeros(3, 2)
    @test b2.store.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.store.velocity_half == zeros(3, 2)
    @test b2.store.acceleration == zeros(3, 2)
    @test b2.store.b_int == zeros(3, 4)
    @test b2.store.b_ext == zeros(3, 2)
    @test b2.store.damage ≈ [2/3, 2/3]
    @test b2.store.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.store.n_active_bonds == [1, 1]

    randpos = rand(3, 4)
    tdh.chunks[2].store.position .= randpos
    Peridynamics.exchange_read_fields!(tdh, 1)

    @test b1.store.position[:,1:2] ≈ [0.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.store.position[:,3:4] ≈ randpos[:,1:2]
    @test b1.store.displacement == zeros(3, 2)
    @test b1.store.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.store.velocity_half == zeros(3, 2)
    @test b1.store.acceleration == zeros(3, 2)
    @test b1.store.b_int == zeros(3, 4)
    @test b1.store.b_ext == zeros(3, 2)
    @test b1.store.damage ≈ [2/3, 2/3]
    @test b1.store.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.store.n_active_bonds == [1, 1]

    @test b2.store.position ≈ randpos
    @test b2.store.displacement == zeros(3, 2)
    @test b2.store.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.store.velocity_half == zeros(3, 2)
    @test b2.store.acceleration == zeros(3, 2)
    @test b2.store.b_int == zeros(3, 4)
    @test b2.store.b_ext == zeros(3, 2)
    @test b2.store.damage ≈ [2/3, 2/3]
    @test b2.store.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.store.n_active_bonds == [1, 1]

    randbint = rand(3, 4)
    b2.store.b_int .= randbint
    b1.store.b_int .+= 1
    Peridynamics.exchange_write_fields!(tdh, 1)

    @test b1.store.position[:,1:2] ≈ [0.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.store.position[:,3:4] ≈ randpos[:,1:2]
    @test b1.store.displacement == zeros(3, 2)
    @test b1.store.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.store.velocity_half == zeros(3, 2)
    @test b1.store.acceleration == zeros(3, 2)
    @test b1.store.b_int[:,1:2] ≈ 1 .+ randbint[:,3:4]
    @test b1.store.b_int[:,3:4] ≈ ones(3, 2)
    @test b1.store.b_ext == zeros(3, 2)
    @test b1.store.damage ≈ [2/3, 2/3]
    @test b1.store.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.store.n_active_bonds == [1, 1]

    @test b2.store.position ≈ randpos
    @test b2.store.displacement == zeros(3, 2)
    @test b2.store.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.store.velocity_half == zeros(3, 2)
    @test b2.store.acceleration == zeros(3, 2)
    @test b2.store.b_int ≈ randbint
    @test b2.store.b_ext == zeros(3, 2)
    @test b2.store.damage ≈ [2/3, 2/3]
    @test b2.store.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.store.n_active_bonds == [1, 1]

    Peridynamics.exchange_write_fields!(tdh, 2)

    @test b1.store.position[:,1:2] ≈ [0.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.store.position[:,3:4] ≈ randpos[:,1:2]
    @test b1.store.displacement == zeros(3, 2)
    @test b1.store.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.store.velocity_half == zeros(3, 2)
    @test b1.store.acceleration == zeros(3, 2)
    @test b1.store.b_int[:,1:2] ≈ 1 .+ randbint[:,3:4]
    @test b1.store.b_int[:,3:4] ≈ ones(3, 2)
    @test b1.store.b_ext == zeros(3, 2)
    @test b1.store.damage ≈ [2/3, 2/3]
    @test b1.store.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.store.n_active_bonds == [1, 1]

    @test b2.store.position ≈ randpos
    @test b2.store.displacement == zeros(3, 2)
    @test b2.store.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.store.velocity_half == zeros(3, 2)
    @test b2.store.acceleration == zeros(3, 2)
    @test b2.store.b_int[:,1:2] ≈ 1 .+ randbint[:,1:2]
    @test b2.store.b_int[:,3:4] ≈ randbint[:,3:4]
    @test b2.store.b_ext == zeros(3, 2)
    @test b2.store.damage ≈ [2/3, 2/3]
    @test b2.store.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.store.n_active_bonds == [1, 1]
end
