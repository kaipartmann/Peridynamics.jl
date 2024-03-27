@testitem "BodyChunk BBMaterial" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    velocity_ic!(body, :a, :x, 1.0)
    velocity_bc!(t->t, body, :a, :x)
    forcedensity_bc!(t->t, body, :a, :x)
    precrack!(body, :a, :b)
    ts = VelocityVerlet(steps=10)
    pd = Peridynamics.PointDecomposition(body, 2)
    bc = Peridynamics.BodyChunk(body, ts, pd, 1)

    @test bc.mat == mat
    @test bc.system isa Peridynamics.BondSystem
    @test bc.system.position == position
    @test bc.system.volume == volume
    @test bc.system.bonds == [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(3, √2, true),
        Peridynamics.Bond(4, √2, true),
    ]
    @test bc.system.n_neighbors == [3, 3]
    @test bc.system.bond_ids == [1:3, 4:6]

    @test bc.ch.point_ids == [1, 2, 3, 4]
    @test bc.ch.loc_points == pd.decomp[1]
    @test bc.ch.halo_points == [3, 4]
    @test bc.ch.hidxs_by_src[2] == 3:4
    for i in 1:4
        @test bc.ch.localizer[i] == i
    end

    #TODO: test the other fields!
end


@testitem "chop_body_threads BBMaterial Val{1}" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    velocity_ic!(body, :a, :x, 1.0)
    velocity_bc!(t->t, body, :a, :x)
    forcedensity_bc!(t->t, body, :a, :x)
    precrack!(body, :a, :b)
    ts = VelocityVerlet(steps=10)
    point_decomp = Peridynamics.PointDecomposition(body, 2)

    body_chunks = Peridynamics.chop_body_threads(body, ts, point_decomp, Val{1}())

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
    @test b1.storage.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.velocity_half == zeros(3, 2)
    @test b1.storage.acceleration == zeros(3, 2)
    @test b1.storage.b_int == zeros(3, 2)
    @test b1.storage.b_ext == zeros(3, 2)
    @test b1.storage.damage ≈ [2/3, 2/3]
    @test b1.storage.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.storage.n_active_bonds == [1, 1]

    @test b1.param isa Peridynamics.BBPointParameters
    @test b1.param.δ ≈ 2.0
    @test b1.param.rho ≈ 1.0
    @test b1.param.E ≈ 1.0
    @test b1.param.nu ≈ 0.25
    @test b1.param.G ≈ 0.4
    @test b1.param.K ≈ 2/3
    @test b1.param.λ ≈ 0.4
    @test b1.param.μ ≈ 0.4
    @test b1.param.Gc ≈ 1.0
    @test b1.param.εc ≈ 0.6454972243679028
    @test b1.param.bc ≈ 0.238732414637843
    @test b1.param == Peridynamics.get_param(b1, 1)
    @test b1.param == Peridynamics.get_param(b1, 2)
    @test b1.param == Peridynamics.get_param(b1, 100)

    @test b1.psets[:a] == [1, 2]
    @test b1.psets[:b] == []

    @test length(b1.sdbcs) == 2
    @test b1.sdbcs[1].fun(0) == 0
    @test b1.sdbcs[1].fun(1) == 1
    @test b1.sdbcs[1].field === :velocity_half
    @test b1.sdbcs[1].point_set === :a
    @test b1.sdbcs[1].dim == 0x01

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

    @test b2.param isa Peridynamics.BBPointParameters
    @test b2.param.δ ≈ 2.0
    @test b2.param.rho ≈ 1.0
    @test b2.param.E ≈ 1.0
    @test b2.param.nu ≈ 0.25
    @test b2.param.G ≈ 0.4
    @test b2.param.K ≈ 2/3
    @test b2.param.λ ≈ 0.4
    @test b2.param.μ ≈ 0.4
    @test b2.param.Gc ≈ 1.0
    @test b2.param.εc ≈ 0.6454972243679028
    @test b2.param.bc ≈ 0.238732414637843
    @test b2.param == Peridynamics.get_param(b2, 1)
    @test b2.param == Peridynamics.get_param(b2, 2)
    @test b2.param == Peridynamics.get_param(b2, 100)

    @test b2.psets[:a] == []
    @test b2.psets[:b] == [1, 2]

    @test length(b2.sdbcs) == 2
    @test b2.sdbcs[1].fun(0) == 0
    @test b2.sdbcs[1].fun(1) == 1
    @test b2.sdbcs[1].field === :velocity_half
    @test b2.sdbcs[1].point_set === :a
    @test b2.sdbcs[1].dim == 0x01

    @test b2.ch.point_ids == [3, 4, 1, 2]
    @test b2.ch.loc_points == [3, 4]
    @test b2.ch.halo_points == [1, 2]
    @test b2.ch.hidxs_by_src[1] == 3:4
    @test b2.ch.localizer[3] == 1
    @test b2.ch.localizer[4] == 2
    @test b2.ch.localizer[1] == 3
    @test b2.ch.localizer[2] == 4
end

@testitem "chop_body_threads BBMaterial Val{N}" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    material!(body, :b, horizon=3, rho=2, E=2, Gc=2)
    velocity_ic!(body, :a, :x, 1.0)
    velocity_bc!(t->t, body, :a, :x)
    forcedensity_bc!(t->t, body, :a, :x)
    precrack!(body, :a, :b)
    ts = VelocityVerlet(steps=10)
    point_decomp = Peridynamics.PointDecomposition(body, 2)

    body_chunks = Peridynamics.chop_body_threads(body, ts, point_decomp, Val{2}())

    b1 = body_chunks[1]
    @test b1 isa Peridynamics.MultiParamBodyChunk
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
    @test b1.storage.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.velocity_half == zeros(3, 2)
    @test b1.storage.acceleration == zeros(3, 2)
    @test b1.storage.b_int == zeros(3, 2)
    @test b1.storage.b_ext == zeros(3, 2)
    @test b1.storage.damage ≈ [2/3, 2/3]
    @test b1.storage.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.storage.n_active_bonds == [1, 1]

    @test length(b1.param) == 2

    pp1 = Peridynamics.get_param(b1, 1)
    @test pp1 isa Peridynamics.BBPointParameters
    @test pp1.δ ≈ 2.0
    @test pp1.rho ≈ 1.0
    @test pp1.E ≈ 1.0
    @test pp1.nu ≈ 0.25
    @test pp1.G ≈ 0.4
    @test pp1.K ≈ 2/3
    @test pp1.λ ≈ 0.4
    @test pp1.μ ≈ 0.4
    @test pp1.Gc ≈ 1.0
    @test pp1.εc ≈ 0.6454972243679028
    @test pp1.bc ≈ 0.238732414637843

    pp2 = Peridynamics.get_param(b1, 2)
    @test pp2 isa Peridynamics.BBPointParameters
    @test pp2.δ ≈ 2.0
    @test pp2.rho ≈ 1.0
    @test pp2.E ≈ 1.0
    @test pp2.nu ≈ 0.25
    @test pp2.G ≈ 0.4
    @test pp2.K ≈ 2/3
    @test pp2.λ ≈ 0.4
    @test pp2.μ ≈ 0.4
    @test pp2.Gc ≈ 1.0
    @test pp2.εc ≈ 0.6454972243679028
    @test pp2.bc ≈ 0.238732414637843

    @test_throws BoundsError Peridynamics.get_param(b1, 3)
    @test_throws BoundsError Peridynamics.get_param(b1, 4)

    @test b1.psets[:a] == [1, 2]
    @test b1.psets[:b] == []

    @test length(b1.sdbcs) == 2
    @test b1.sdbcs[1].fun(0) == 0
    @test b1.sdbcs[1].fun(1) == 1
    @test b1.sdbcs[1].field === :velocity_half
    @test b1.sdbcs[1].point_set === :a
    @test b1.sdbcs[1].dim == 0x01

    @test b1.ch.point_ids == [1, 2, 3, 4]
    @test b1.ch.loc_points == 1:2
    @test b1.ch.halo_points == [3, 4]
    @test b1.ch.hidxs_by_src[2] == 3:4
    for i in 1:4
        @test b1.ch.localizer[i] == i
    end

    b2 = body_chunks[2]
    @test b2 isa Peridynamics.MultiParamBodyChunk
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

    @test length(b2.param) == 2

    pp3 = Peridynamics.get_param(b2, 1)
    @test pp3 isa Peridynamics.BBPointParameters
    @test pp3.δ ≈ 3.0
    @test pp3.rho ≈ 2.0
    @test pp3.E ≈ 2.0
    @test pp3.nu ≈ 0.25
    @test pp3.G ≈ 0.8
    @test pp3.K ≈ 4/3
    @test pp3.λ ≈ 0.8
    @test pp3.μ ≈ 0.8
    @test pp3.Gc ≈ 2.0
    @test pp3.εc ≈ 0.5270462766947299
    @test pp3.bc ≈ 0.0943140403507528

    pp4= Peridynamics.get_param(b2, 2)
    @test pp4 isa Peridynamics.BBPointParameters
    @test pp4.δ ≈ 3.0
    @test pp4.rho ≈ 2.0
    @test pp4.E ≈ 2.0
    @test pp4.nu ≈ 0.25
    @test pp4.G ≈ 0.8
    @test pp4.K ≈ 4/3
    @test pp4.λ ≈ 0.8
    @test pp4.μ ≈ 0.8
    @test pp4.Gc ≈ 2.0
    @test pp4.εc ≈ 0.5270462766947299
    @test pp4.bc ≈ 0.0943140403507528

    @test_throws BoundsError Peridynamics.get_param(b1, 3)
    @test_throws BoundsError Peridynamics.get_param(b1, 4)

    @test b2.psets[:a] == []
    @test b2.psets[:b] == [1, 2]

    @test length(b2.sdbcs) == 2
    @test b2.sdbcs[1].fun(0) == 0
    @test b2.sdbcs[1].fun(1) == 1
    @test b2.sdbcs[1].field === :velocity_half
    @test b2.sdbcs[1].point_set === :a
    @test b2.sdbcs[1].dim == 0x01

    @test b2.ch.point_ids == [3, 4, 1, 2]
    @test b2.ch.loc_points == [3, 4]
    @test b2.ch.halo_points == [1, 2]
    @test b2.ch.hidxs_by_src[1] == 3:4
    @test b2.ch.localizer[3] == 1
    @test b2.ch.localizer[4] == 2
    @test b2.ch.localizer[1] == 3
    @test b2.ch.localizer[2] == 4
end
