@testitem "precrack bond chunk" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    precrack!(body, :a, :b)
    ts = VelocityVerlet(steps=10)
    point_decomp = Peridynamics.PointDecomposition(body, 2)

    body_chunks = Peridynamics.chop_body_threads(body, ts, point_decomp)

    b1 = body_chunks[1]
    @test b1 isa Peridynamics.BodyChunk
    @test b1.mat == BBMaterial()
    @test b1.system isa Peridynamics.BondSystem
    @test b1.system.position == position
    @test b1.system.volume == volume
    @test b1.system.bonds == [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(3, √2, true),
        Peridynamics.Bond(4, √2, true),
    ]
    @test b1.system.n_neighbors == [3, 3]
    @test b1.system.bond_ids == [1:3, 4:6]

    @test b1.storage.position == position
    @test b1.storage.displacement == zeros(3, 2)
    @test b1.storage.velocity == zeros(3, 2)
    @test b1.storage.velocity_half == zeros(3, 2)
    @test b1.storage.acceleration == zeros(3, 2)
    @test b1.storage.b_int == zeros(3, 2)
    @test b1.storage.b_ext == zeros(3, 2)
    @test b1.storage.damage ≈ [2/3, 2/3]
    @test b1.storage.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.storage.n_active_bonds == [1, 1]

    @test b1.ch.point_ids == [1, 2, 3, 4]
    @test b1.ch.loc_points == 1:2
    @test b1.ch.halo_points == [3, 4]
    @test b1.ch.hidxs_by_src[2] == 3:4
    for i in 1:4
        @test b1.ch.localizer[i] == i
    end

    b2 = body_chunks[2]
    @test b2 isa Peridynamics.BodyChunk
    @test b2.mat == BBMaterial()
    @test b2.system isa Peridynamics.BondSystem
    @test b2.system.position == position[:, [3, 4, 1, 2]]
    @test b2.system.volume == volume[[3, 4, 1, 2]]
    @test b2.system.bonds == [
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, √2, true),
        Peridynamics.Bond(2, √2, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, √2, true),
        Peridynamics.Bond(1, √2, true),
    ]
    @test b2.system.n_neighbors == [3, 3]
    @test b2.system.bond_ids == [1:3, 4:6]

    @test b2.storage.position == position[:, [3, 4, 1, 2]]
    @test b2.storage.displacement == zeros(3, 2)
    @test b2.storage.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.storage.velocity_half == zeros(3, 2)
    @test b2.storage.acceleration == zeros(3, 2)
    @test b2.storage.b_int == zeros(3, 2)
    @test b2.storage.b_ext == zeros(3, 2)
    @test b2.storage.damage ≈ [2/3, 2/3]
    @test b2.storage.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.storage.n_active_bonds == [1, 1]

    @test b2.ch.point_ids == [3, 4, 1, 2]
    @test b2.ch.loc_points == [3, 4]
    @test b2.ch.halo_points == [1, 2]
    @test b2.ch.hidxs_by_src[1] == 3:4
    @test b2.ch.localizer[3] == 1
    @test b2.ch.localizer[4] == 2
    @test b2.ch.localizer[1] == 3
    @test b2.ch.localizer[2] == 4
end

@testitem "precrack bond chunk with bond filter" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    precrack!(body, :a, :b; update_dmg=false)
    ts = VelocityVerlet(steps=10)
    point_decomp = Peridynamics.PointDecomposition(body, 2)

    body_chunks = Peridynamics.chop_body_threads(body, ts, point_decomp)

    b1 = body_chunks[1]
    @test b1 isa Peridynamics.BodyChunk
    @test b1.mat == BBMaterial()
    @test b1.system isa Peridynamics.BondSystem
    @test b1.system.position == position[:, 1:2]
    @test b1.system.volume == volume[1:2]
    @test b1.system.bonds == [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
    ]
    @test b1.system.n_neighbors == [1, 1]
    @test b1.system.bond_ids == [1:1, 2:2]

    @test b1.storage.position == position[:, 1:2]
    @test b1.storage.displacement == zeros(3, 2)
    @test b1.storage.velocity == zeros(3, 2)
    @test b1.storage.velocity_half == zeros(3, 2)
    @test b1.storage.acceleration == zeros(3, 2)
    @test b1.storage.b_int == zeros(3, 2)
    @test b1.storage.b_ext == zeros(3, 2)
    @test b1.storage.damage == zeros(2)
    @test b1.storage.bond_active == [1, 1]
    @test b1.storage.n_active_bonds == [1, 1]

    @test b1.ch.point_ids == [1, 2]
    @test b1.ch.loc_points == 1:2
    @test b1.ch.halo_points == Int[]
    @test b1.ch.hidxs_by_src == Dict{Int,UnitRange{Int}}()
    for i in 1:2
        @test b1.ch.localizer[i] == i
    end

    b2 = body_chunks[2]
    @test b2 isa Peridynamics.BodyChunk
    @test b2.mat == BBMaterial()
    @test b2.system isa Peridynamics.BondSystem
    @test b2.system.position == position[:, 3:4]
    @test b2.system.volume == volume[3:4]
    @test b2.system.bonds == [
        Peridynamics.Bond(2, √2, true),
        Peridynamics.Bond(1, √2, true),
    ]
    @test b2.system.n_neighbors == [1, 1]
    @test b2.system.bond_ids == [1:1, 2:2]

    @test b2.storage.position == position[:, 3:4]
    @test b2.storage.displacement == zeros(3, 2)
    @test b2.storage.velocity == zeros(3, 2)
    @test b2.storage.velocity_half == zeros(3, 2)
    @test b2.storage.acceleration == zeros(3, 2)
    @test b2.storage.b_int == zeros(3, 2)
    @test b2.storage.b_ext == zeros(3, 2)
    @test b2.storage.damage == zeros(2)
    @test b2.storage.bond_active == [1, 1]
    @test b2.storage.n_active_bonds == [1, 1]

    @test b2.ch.point_ids == [3, 4]
    @test b2.ch.loc_points == [3, 4]
    @test b2.ch.halo_points == Int[]
    @test b2.ch.hidxs_by_src == Dict{Int,UnitRange{Int}}()
    @test b2.ch.localizer[3] == 1
    @test b2.ch.localizer[4] == 2
end

@testitem "precrack bond system with filter" begin
    # setup
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    point_set!(body, :set_a, 1:2)
    point_set!(body, :set_b, 3:4)
    precrack!(body, :set_a, :set_b; update_dmg=false)
    pd = Peridynamics.PointDecomposition(body, 2)

    # 1
    bs, ch = Peridynamics.BondSystem(body, pd, 1)

    @test bs.position == position[:, 1:2]
    @test bs.volume == volume[1:2]
    @test bs.bonds == [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
    ]
    @test bs.n_neighbors == [1, 1]
    @test bs.bond_ids == [1:1, 2:2]

    @test ch.point_ids == [1, 2]
    @test ch.loc_points == [1, 2]
    @test ch.halo_points == Int[]
    @test ch.hidxs_by_src == Dict{Int,UnitRange{Int}}()

    for i in 1:2
        @test ch.localizer[i] == i
    end

    # 2
    bs, ch = Peridynamics.BondSystem(body, pd, 2)

    @test bs.position == position[:, 3:4]
    @test bs.volume == volume[3:4]
    @test bs.bonds == [
        Peridynamics.Bond(2, √2, true),
        Peridynamics.Bond(1, √2, true),
    ]
    @test bs.n_neighbors == [1, 1]
    @test bs.bond_ids == [1:1, 2:2]

    @test ch.point_ids == [3, 4]
    @test ch.loc_points == [3, 4]
    @test ch.halo_points == Int[]
    @test ch.hidxs_by_src == Dict{Int,UnitRange{Int}}()
    @test ch.localizer[3] == 1
    @test ch.localizer[4] == 2
end
