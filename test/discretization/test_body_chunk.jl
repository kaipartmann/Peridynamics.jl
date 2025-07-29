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
    ps = Peridynamics.get_param_spec(body)
    chunk = Peridynamics.BodyChunk(body, ts, pd, 1, ps)

    @test chunk.mat == mat
    @test chunk.system isa Peridynamics.BondSystem
    @test chunk.system.position == position
    @test chunk.system.volume == volume
    @test chunk.system.bonds == [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(3, √2, true),
        Peridynamics.Bond(4, √2, true),
    ]
    @test chunk.system.n_neighbors == [3, 3]
    @test chunk.system.bond_ids == [1:3, 4:6]

    ch = chunk.system.chunk_handler
    @test ch.point_ids == [1, 2, 3, 4]
    @test ch.loc_points == pd.decomp[1]
    @test ch.halo_points == [3, 4]
    @test ch.hidxs_by_src[2] == 3:4
    for i in 1:4
        @test ch.localizer[i] == i
    end
    @test Peridynamics.n_loc_points(chunk) == 2
    @test Peridynamics.n_points(chunk) == 4

    #TODO: test the other fields!
end


@testitem "chop_body_threads BBMaterial N=1" begin
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
    param_spec = Peridynamics.get_param_spec(body)

    body_chunks = Peridynamics.chop_body_threads(body, ts, point_decomp, param_spec)

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

    @test b1.paramsetup isa Peridynamics.StandardPointParameters
    b1_params = Peridynamics.get_params(b1, 1)
    @test b1_params.δ ≈ 2.0
    @test b1_params.rho ≈ 1.0
    @test b1_params.E ≈ 1.0
    @test b1_params.nu ≈ 0.25
    @test b1_params.G ≈ 0.4
    @test b1_params.K ≈ 2/3
    @test b1_params.λ ≈ 0.4
    @test b1_params.μ ≈ 0.4
    @test b1_params.Gc ≈ 1.0
    @test b1_params.εc ≈ 0.6454972243679028
    @test b1_params.bc ≈ 0.238732414637843 * π / (8 * 0.9605919564548167) # correction factor for cube neighborhood
    @test b1_params == Peridynamics.get_params(b1, 2)
    @test b1_params == Peridynamics.get_params(b1, 100)

    @test b1.psets[:a] == [1, 2]
    @test b1.psets[:b] == []

    @test length(b1.sdbcs) == 2
    @test b1.sdbcs[1].fun(0) == 0
    @test b1.sdbcs[1].fun(1) == 1
    @test b1.sdbcs[1].field === :velocity_half
    @test b1.sdbcs[1].point_set === :a
    @test b1.sdbcs[1].dim == 0x01

    ch1 = b1.system.chunk_handler
    @test ch1.point_ids == [1, 2, 3, 4]
    @test ch1.loc_points == 1:2
    @test ch1.halo_points == [3, 4]
    @test ch1.hidxs_by_src[2] == 3:4
    for i in 1:4
        @test ch1.localizer[i] == i
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

    @test b2.paramsetup isa Peridynamics.StandardPointParameters
    b2_params = Peridynamics.get_params(b2, 1)
    @test b2_params.δ ≈ 2.0
    @test b2_params.rho ≈ 1.0
    @test b2_params.E ≈ 1.0
    @test b2_params.nu ≈ 0.25
    @test b2_params.G ≈ 0.4
    @test b2_params.K ≈ 2/3
    @test b2_params.λ ≈ 0.4
    @test b2_params.μ ≈ 0.4
    @test b2_params.Gc ≈ 1.0
    @test b2_params.εc ≈ 0.6454972243679028
    @test b2_params.bc ≈ 0.238732414637843 * π / (8 * 0.9605919564548167) # correction factor for cube neighborhood
    @test b2_params == Peridynamics.get_params(b2, 2)
    @test b2_params == Peridynamics.get_params(b2, 100)

    @test b2.psets[:a] == []
    @test b2.psets[:b] == [1, 2]

    @test length(b2.sdbcs) == 2
    @test b2.sdbcs[1].fun(0) == 0
    @test b2.sdbcs[1].fun(1) == 1
    @test b2.sdbcs[1].field === :velocity_half
    @test b2.sdbcs[1].point_set === :a
    @test b2.sdbcs[1].dim == 0x01

    ch2 = b2.system.chunk_handler
    @test ch2.point_ids == [3, 4, 1, 2]
    @test ch2.loc_points == [3, 4]
    @test ch2.halo_points == [1, 2]
    @test ch2.hidxs_by_src[1] == 3:4
    @test ch2.localizer[3] == 1
    @test ch2.localizer[4] == 2
    @test ch2.localizer[1] == 3
    @test ch2.localizer[2] == 4
end

@testitem "chop_body_threads BBMaterial N=2" begin
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
    param_spec = Peridynamics.get_param_spec(body)

    body_chunks = Peridynamics.chop_body_threads(body, ts, point_decomp, param_spec)

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

    @test length(b1.paramsetup.parameters) == 2

    pp1 = Peridynamics.get_params(b1, 1)
    @test pp1 isa Peridynamics.StandardPointParameters
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
    @test pp1.bc ≈ 0.238732414637843 * π / (8 * 0.9605919564548167) # correction factor for cube neighborhood

    pp2 = Peridynamics.get_params(b1, 2)
    @test pp2 isa Peridynamics.StandardPointParameters
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
    @test pp2.bc ≈ 0.238732414637843 * π / (8 * 0.9605919564548167) # correction factor for cube neighborhood

    pp3 = Peridynamics.get_params(b1, 3)
    @test pp3 isa Peridynamics.StandardPointParameters
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
    @test pp3.bc ≈ 0.0943140403507528 * π / (8 * 0.9605919564548167) # correction factor for cube neighborhood

    pp4= Peridynamics.get_params(b1, 4)
    @test pp4 isa Peridynamics.StandardPointParameters
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
    @test pp4.bc ≈ 0.0943140403507528 * π / (8 * 0.9605919564548167) # correction factor for cube neighborhood

    @test b1.psets[:a] == [1, 2]
    @test b1.psets[:b] == []

    @test length(b1.sdbcs) == 2
    @test b1.sdbcs[1].fun(0) == 0
    @test b1.sdbcs[1].fun(1) == 1
    @test b1.sdbcs[1].field === :velocity_half
    @test b1.sdbcs[1].point_set === :a
    @test b1.sdbcs[1].dim == 0x01

    ch1 = b1.system.chunk_handler
    @test ch1.point_ids == [1, 2, 3, 4]
    @test ch1.loc_points == 1:2
    @test ch1.halo_points == [3, 4]
    @test ch1.hidxs_by_src[2] == 3:4
    for i in 1:4
        @test ch1.localizer[i] == i
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

    @test length(b2.paramsetup.parameters) == 2

    pp3 = Peridynamics.get_params(b2, 1)
    @test pp3 isa Peridynamics.StandardPointParameters
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
    @test pp3.bc ≈ 0.0943140403507528 * π / (8 * 0.9605919564548167) # correction factor for cube neighborhood

    pp4 = Peridynamics.get_params(b2, 2)
    @test pp4 isa Peridynamics.StandardPointParameters
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
    @test pp4.bc ≈ 0.0943140403507528 * π / (8 * 0.9605919564548167) # correction factor for cube neighborhood

    pp1 = Peridynamics.get_params(b2, 3)
    @test pp1 isa Peridynamics.StandardPointParameters
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
    @test pp1.bc ≈ 0.238732414637843 * π / (8 * 0.9605919564548167) # correction factor for cube neighborhood

    pp2 = Peridynamics.get_params(b2, 4)
    @test pp2 isa Peridynamics.StandardPointParameters
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
    @test pp2.bc ≈ 0.238732414637843 * π / (8 * 0.9605919564548167) # correction factor for cube neighborhood

    @test b2.psets[:a] == []
    @test b2.psets[:b] == [1, 2]

    @test length(b2.sdbcs) == 2
    @test b2.sdbcs[1].fun(0) == 0
    @test b2.sdbcs[1].fun(1) == 1
    @test b2.sdbcs[1].field === :velocity_half
    @test b2.sdbcs[1].point_set === :a
    @test b2.sdbcs[1].dim == 0x01

    ch2 = b2.system.chunk_handler
    @test ch2.point_ids == [3, 4, 1, 2]
    @test ch2.loc_points == [3, 4]
    @test ch2.halo_points == [1, 2]
    @test ch2.hidxs_by_src[1] == 3:4
    @test ch2.localizer[3] == 1
    @test ch2.localizer[4] == 2
    @test ch2.localizer[1] == 3
    @test ch2.localizer[2] == 4
end
