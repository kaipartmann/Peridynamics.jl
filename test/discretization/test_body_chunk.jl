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
    @test Peridynamics.get_n_loc_points(chunk) == 2
    @test Peridynamics.get_n_points(chunk) == 4

    #TODO: test the other fields!
end

@testitem "BodyChunk, 10 points, 2 chunks" begin
    position = [0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    volume = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=1.1, rho=1, E=1)
    ts = VelocityVerlet(steps=10)
    pd = Peridynamics.PointDecomposition(body, 2)
    ps = Peridynamics.get_param_spec(body)

    # first chunk
    chunk = Peridynamics.BodyChunk(body, ts, pd, 1, ps)

    @test chunk.mat == mat
    @test chunk.system isa Peridynamics.BondSystem
    point_ids1 = [1, 2, 3, 4, 5, 6]
    @test chunk.system.position == position[:, point_ids1]
    @test chunk.system.volume == volume[point_ids1]
    @test chunk.system.bonds == [
        Peridynamics.Bond(2, 1.0, false),
        Peridynamics.Bond(1, 1.0, false),
        Peridynamics.Bond(3, 1.0, false),
        Peridynamics.Bond(2, 1.0, false),
        Peridynamics.Bond(4, 1.0, false),
        Peridynamics.Bond(3, 1.0, false),
        Peridynamics.Bond(5, 1.0, false),
        Peridynamics.Bond(4, 1.0, false),
        Peridynamics.Bond(6, 1.0, false),
    ]
    @test chunk.system.n_neighbors == [1, 2, 2, 2, 2]
    @test chunk.system.bond_ids == [1:1, 2:3, 4:5, 6:7, 8:9]

    ch = chunk.system.chunk_handler
    @test ch.n_loc_points == 5
    @test ch.point_ids == [1, 2, 3, 4, 5, 6]
    @test ch.loc_points == 1:5
    @test ch.halo_points == [6]
    @test keys(ch.hidxs_by_src) == Set([2])
    @test ch.hidxs_by_src[2] == 6:6
    @test keys(ch.localizer) == Set(1:6)
    for i in 1:6
        @test ch.localizer[i] == i
    end

    @test size(chunk.storage.displacement) == (3, 5)
    @test Peridynamics.get_loc_point_data(chunk.storage, chunk.system, :displacement) ≈
          chunk.storage.displacement[:, 1:5]

    @test Peridynamics.each_point_idx(chunk.system) == 1:5
    @test collect(Peridynamics.each_loc_dof(chunk))[:] == [1:3:13; 2:3:14; 3:3:15]
    @test collect(Peridynamics.each_dof(chunk))[:] == [1:3:16; 2:3:17; 3:3:18]
    @test setdiff(chunk.condhandler.free_dofs,
                  collect(Peridynamics.each_dof(chunk))[:]) == Int[]
    @test Peridynamics.get_n_loc_dof(chunk) == 15
    @test Peridynamics.get_n_dof(chunk) == 18

    # second chunk
    chunk = Peridynamics.BodyChunk(body, ts, pd, 2, ps)

    @test chunk.mat == mat
    @test chunk.system isa Peridynamics.BondSystem
    point_ids2 = [6, 7, 8, 9, 10, 5]
    @test chunk.system.position == position[:, point_ids2]
    @test chunk.system.volume == volume[point_ids2]
    @test chunk.system.bonds == [
        Peridynamics.Bond(6, 1.0, false),
        Peridynamics.Bond(2, 1.0, false),
        Peridynamics.Bond(1, 1.0, false),
        Peridynamics.Bond(3, 1.0, false),
        Peridynamics.Bond(2, 1.0, false),
        Peridynamics.Bond(4, 1.0, false),
        Peridynamics.Bond(3, 1.0, false),
        Peridynamics.Bond(5, 1.0, false),
        Peridynamics.Bond(4, 1.0, false),
    ]
    @test chunk.system.n_neighbors == [2, 2, 2, 2, 1]
    @test chunk.system.bond_ids == [1:2, 3:4, 5:6, 7:8, 9:9]

    ch = chunk.system.chunk_handler
    @test ch.n_loc_points == 5
    @test ch.point_ids == [6, 7, 8, 9, 10, 5]
    @test ch.loc_points == 6:10
    @test ch.halo_points == [5]
    @test keys(ch.hidxs_by_src) == Set([1])
    @test ch.hidxs_by_src[1] == 6:6
    @test keys(ch.localizer) == Set(5:10)
    @test ch.localizer[6] == 1
    @test ch.localizer[7] == 2
    @test ch.localizer[8] == 3
    @test ch.localizer[9] == 4
    @test ch.localizer[10] == 5
    @test ch.localizer[5] == 6

    @test size(chunk.storage.displacement) == (3, 5) # 6 sec loc points
    @test Peridynamics.get_loc_point_data(chunk.storage, chunk.system, :displacement) ≈
          chunk.storage.displacement[:, 1:5]

    @test Peridynamics.each_point_idx(chunk.system) == 1:5
    @test collect(Peridynamics.each_loc_dof(chunk.system))[:] == [1:3:13; 2:3:14; 3:3:15]
    @test collect(Peridynamics.each_dof(chunk.system))[:] == [1:3:16; 2:3:17; 3:3:18]
    @test setdiff(chunk.condhandler.free_dofs,
                  collect(Peridynamics.each_dof(chunk.system))[:]) == Int[]
end

@testitem "BodyChunk, 10 points, 5 chunks" begin
    position = [0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    volume = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=1.01, rho=1, E=1)
    ts = VelocityVerlet(steps=10)
    pd = Peridynamics.PointDecomposition(body, 5)
    ps = Peridynamics.get_param_spec(body)

    # chunk 1
    chunk = Peridynamics.BodyChunk(body, ts, pd, 1, ps)

    @test chunk.mat == mat
    @test chunk.system isa Peridynamics.BondSystem
    point_ids1 = [1, 2, 3]
    @test chunk.system.position == position[:, point_ids1]
    @test chunk.system.volume == volume[point_ids1]
    @test chunk.system.bonds == [
        Peridynamics.Bond(2, 1.0, false),
        Peridynamics.Bond(1, 1.0, false),
        Peridynamics.Bond(3, 1.0, false),
    ]
    @test chunk.system.n_neighbors == [1, 2]
    @test chunk.system.bond_ids == [1:1, 2:3]

    ch = chunk.system.chunk_handler
    @test ch.n_loc_points == 2
    @test ch.point_ids == [1, 2, 3]
    @test ch.loc_points == 1:2
    @test ch.halo_points == [3]
    @test keys(ch.hidxs_by_src) == Set([2])
    @test ch.hidxs_by_src[2] == 3:3
    @test keys(ch.localizer) == Set(1:3)
    for i in 1:3
        @test ch.localizer[i] == i
    end

    @test size(chunk.storage.displacement) == (3, 2) # 3 sec loc points
    @test Peridynamics.get_loc_point_data(chunk.storage, chunk.system, :displacement) ≈
          chunk.storage.displacement[:, 1:2]

    @test Peridynamics.each_point_idx(chunk.system) == 1:2
    @test collect(Peridynamics.each_loc_dof(chunk.system))[:] == [1, 4, 2, 5, 3, 6]
    @test collect(Peridynamics.each_dof(chunk.system))[:] == [1:3:7; 2:3:8; 3:3:9]
    @test setdiff(chunk.condhandler.free_dofs,
                  collect(Peridynamics.each_dof(chunk.system))[:]) == Int[]

    # chunk 2
    chunk = Peridynamics.BodyChunk(body, ts, pd, 2, ps)

    @test chunk.mat == mat
    @test chunk.system isa Peridynamics.BondSystem
    point_ids2 = [3, 4, 2, 5]
    @test chunk.system.position == position[:, point_ids2]
    @test chunk.system.volume == volume[point_ids2]
    @test chunk.system.bonds == [
        Peridynamics.Bond(3, 1.0, false),
        Peridynamics.Bond(2, 1.0, false),
        Peridynamics.Bond(1, 1.0, false),
        Peridynamics.Bond(4, 1.0, false),
    ]
    @test chunk.system.n_neighbors == [2, 2]
    @test chunk.system.bond_ids == [1:2, 3:4]

    ch = chunk.system.chunk_handler
    @test ch.n_loc_points == 2
    @test ch.point_ids == [3, 4, 2, 5]
    @test ch.loc_points == 3:4
    @test ch.halo_points == [2, 5]
    @test keys(ch.hidxs_by_src) == Set([1, 3])
    @test ch.hidxs_by_src[1] == 3:3
    @test ch.hidxs_by_src[3] == 4:4
    @test keys(ch.localizer) == Set(2:5)
    @test ch.localizer[3] == 1
    @test ch.localizer[4] == 2
    @test ch.localizer[2] == 3
    @test ch.localizer[5] == 4

    @test size(chunk.storage.displacement) == (3, 2) # 4 sec loc points
    @test Peridynamics.get_loc_point_data(chunk.storage, chunk.system, :displacement) ≈
          chunk.storage.displacement[:, 1:2]

    @test Peridynamics.each_point_idx(chunk.system) == 1:2
    @test collect(Peridynamics.each_loc_dof(chunk.system))[:] == [1, 4, 2, 5, 3, 6]
    @test collect(Peridynamics.each_dof(chunk.system))[:] == [1:3:10; 2:3:11; 3:3:12]
    @test setdiff(chunk.condhandler.free_dofs,
                  collect(Peridynamics.each_dof(chunk.system))[:]) == Int[]

    # chunk 3
    chunk = Peridynamics.BodyChunk(body, ts, pd, 3, ps)

    @test chunk.mat == mat
    @test chunk.system isa Peridynamics.BondSystem
    point_ids3 = [5, 6, 4, 7]
    @test chunk.system.position == position[:, point_ids3]
    @test chunk.system.volume == volume[point_ids3]
    @test chunk.system.bonds == [
        Peridynamics.Bond(3, 1.0, false),
        Peridynamics.Bond(2, 1.0, false),
        Peridynamics.Bond(1, 1.0, false),
        Peridynamics.Bond(4, 1.0, false),
    ]
    @test chunk.system.n_neighbors == [2, 2]
    @test chunk.system.bond_ids == [1:2, 3:4]

    ch = chunk.system.chunk_handler
    @test ch.n_loc_points == 2
    @test ch.point_ids == [5, 6, 4, 7]
    @test ch.loc_points == 5:6
    @test ch.halo_points == [4, 7]
    @test keys(ch.hidxs_by_src) == Set([2, 4])
    @test ch.hidxs_by_src[2] == 3:3
    @test ch.hidxs_by_src[4] == 4:4
    @test keys(ch.localizer) == Set(4:7)
    @test ch.localizer[5] == 1
    @test ch.localizer[6] == 2
    @test ch.localizer[4] == 3
    @test ch.localizer[7] == 4

    @test size(chunk.storage.displacement) == (3, 2) # 4 sec loc points
    @test Peridynamics.get_loc_point_data(chunk.storage, chunk.system, :displacement) ≈
          chunk.storage.displacement[:, 1:2]

    @test Peridynamics.each_point_idx(chunk.system) == 1:2
    @test collect(Peridynamics.each_loc_dof(chunk.system))[:] == [1, 4, 2, 5, 3, 6]
    @test collect(Peridynamics.each_dof(chunk.system))[:] == [1:3:10; 2:3:11; 3:3:12]
    @test setdiff(chunk.condhandler.free_dofs,
                  collect(Peridynamics.each_dof(chunk.system))[:]) == Int[]

    # chunk 4
    chunk = Peridynamics.BodyChunk(body, ts, pd, 4, ps)

    @test chunk.mat == mat
    @test chunk.system isa Peridynamics.BondSystem
    point_ids4 = [7, 8, 6, 9]
    @test chunk.system.position == position[:, point_ids4]
    @test chunk.system.volume == volume[point_ids4]
    @test chunk.system.bonds == [
        Peridynamics.Bond(3, 1.0, false),
        Peridynamics.Bond(2, 1.0, false),
        Peridynamics.Bond(1, 1.0, false),
        Peridynamics.Bond(4, 1.0, false),
    ]
    @test chunk.system.n_neighbors == [2, 2]
    @test chunk.system.bond_ids == [1:2, 3:4]

    ch = chunk.system.chunk_handler
    @test ch.n_loc_points == 2
    @test ch.point_ids == [7, 8, 6, 9]
    @test ch.loc_points == 7:8
    @test ch.halo_points == [6, 9]
    @test keys(ch.hidxs_by_src) == Set([3, 5])
    @test ch.hidxs_by_src[3] == 3:3
    @test ch.hidxs_by_src[5] == 4:4
    @test keys(ch.localizer) == Set(6:9)
    @test ch.localizer[7] == 1
    @test ch.localizer[8] == 2
    @test ch.localizer[6] == 3
    @test ch.localizer[9] == 4

    @test size(chunk.storage.displacement) == (3, 2) # 4 sec loc points
    @test Peridynamics.get_loc_point_data(chunk.storage, chunk.system, :displacement) ≈
          chunk.storage.displacement[:, 1:2]

    @test Peridynamics.each_point_idx(chunk.system) == 1:2
    @test collect(Peridynamics.each_loc_dof(chunk.system))[:] == [1, 4, 2, 5, 3, 6]
    @test collect(Peridynamics.each_dof(chunk.system))[:] == [1:3:10; 2:3:11; 3:3:12]
    @test setdiff(chunk.condhandler.free_dofs,
                  collect(Peridynamics.each_dof(chunk.system))[:]) == Int[]

    # chunk 5
    chunk = Peridynamics.BodyChunk(body, ts, pd, 5, ps)

    @test chunk.mat == mat
    @test chunk.system isa Peridynamics.BondSystem
    point_ids5 = [9, 10, 8]
    @test chunk.system.position == position[:, point_ids5]
    @test chunk.system.volume == volume[point_ids5]
    @test chunk.system.bonds == [
        Peridynamics.Bond(3, 1.0, false),
        Peridynamics.Bond(2, 1.0, false),
        Peridynamics.Bond(1, 1.0, false),
    ]
    @test chunk.system.n_neighbors == [2, 1]
    @test chunk.system.bond_ids == [1:2, 3:3]

    ch = chunk.system.chunk_handler
    @test ch.n_loc_points == 2
    @test ch.point_ids == [9, 10, 8]
    @test ch.loc_points == 9:10
    @test ch.halo_points == [8]
    @test keys(ch.hidxs_by_src) == Set([4])
    @test ch.hidxs_by_src[4] == 3:3
    @test keys(ch.localizer) == Set(8:10)
    @test ch.localizer[9] == 1
    @test ch.localizer[10] == 2
    @test ch.localizer[8] == 3

    @test size(chunk.storage.displacement) == (3, 2) # 3 sec loc points
    @test Peridynamics.get_loc_point_data(chunk.storage, chunk.system, :displacement) ≈
          chunk.storage.displacement[:, 1:2]

    @test Peridynamics.each_point_idx(chunk.system) == 1:2
    @test collect(Peridynamics.each_loc_dof(chunk.system))[:] == [1, 4, 2, 5, 3, 6]
    @test collect(Peridynamics.each_dof(chunk.system))[:] == [1:3:7; 2:3:8; 3:3:9]
    @test setdiff(chunk.condhandler.free_dofs,
                  collect(Peridynamics.each_dof(chunk.system))[:]) == Int[]
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
    @test b1_params.bc ≈ 0.238732414637843
    @test b1_params == Peridynamics.get_params(b1, 2)
    @test b1_params == Peridynamics.get_params(b1, 100)

    condhandler1 = b1.condhandler
    @test condhandler1.loc_point_sets[:a] == [1, 2]
    @test condhandler1.loc_point_sets[:b] == []
    @test length(condhandler1.single_dim_bcs) == 2
    @test condhandler1.single_dim_bcs[1].fun(0) == 0
    @test condhandler1.single_dim_bcs[1].fun(1) == 1
    @test condhandler1.single_dim_bcs[1].field === :velocity_half
    @test condhandler1.single_dim_bcs[1].point_set === :a
    @test condhandler1.single_dim_bcs[1].dim == 0x01

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
    @test b2_params.bc ≈ 0.238732414637843
    @test b2_params == Peridynamics.get_params(b2, 2)
    @test b2_params == Peridynamics.get_params(b2, 100)

    condhandler2 = b2.condhandler
    @test condhandler2.loc_point_sets[:a] == []
    @test condhandler2.loc_point_sets[:b] == [1, 2]
    @test length(condhandler2.single_dim_bcs) == 2
    @test condhandler2.single_dim_bcs[1].fun(0) == 0
    @test condhandler2.single_dim_bcs[1].fun(1) == 1
    @test condhandler2.single_dim_bcs[1].field === :velocity_half
    @test condhandler2.single_dim_bcs[1].point_set === :a
    @test condhandler2.single_dim_bcs[1].dim == 0x01

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
    @test pp1.bc ≈ 0.238732414637843

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
    @test pp2.bc ≈ 0.238732414637843

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
    @test pp3.bc ≈ 0.0943140403507528

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
    @test pp4.bc ≈ 0.0943140403507528

    condhandler1 = b1.condhandler
    @test condhandler1.loc_point_sets[:a] == [1, 2]
    @test condhandler1.loc_point_sets[:b] == []
    @test length(condhandler1.single_dim_bcs) == 2
    @test condhandler1.single_dim_bcs[1].fun(0) == 0
    @test condhandler1.single_dim_bcs[1].fun(1) == 1
    @test condhandler1.single_dim_bcs[1].field === :velocity_half
    @test condhandler1.single_dim_bcs[1].point_set === :a
    @test condhandler1.single_dim_bcs[1].dim == 0x01

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
    @test pp3.bc ≈ 0.0943140403507528

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
    @test pp4.bc ≈ 0.0943140403507528

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
    @test pp1.bc ≈ 0.238732414637843

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
    @test pp2.bc ≈ 0.238732414637843

    condhandler2 = b2.condhandler
    @test condhandler2.loc_point_sets[:a] == []
    @test condhandler2.loc_point_sets[:b] == [1, 2]
    @test length(condhandler2.single_dim_bcs) == 2
    @test condhandler2.single_dim_bcs[1].fun(0) == 0
    @test condhandler2.single_dim_bcs[1].fun(1) == 1
    @test condhandler2.single_dim_bcs[1].field === :velocity_half
    @test condhandler2.single_dim_bcs[1].point_set === :a
    @test condhandler2.single_dim_bcs[1].dim == 0x01

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

@testitem "chop_body_threads BBMaterial VelocityVerlet 10-points N=2" begin
    position = [0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    volume = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=1.01, rho=1, E=1, Gc=1)
    point_set!(body, :a, 1:5)
    point_set!(body, :b, 6:10)
    material!(body, :b, horizon=2.01, rho=2, E=2, Gc=2)
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
    @test b1.system.position == position[:, 1:6]
    @test b1.system.volume == volume[1:6]
    @test b1.system.bonds == [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(4, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(5, 1.0, true),
        Peridynamics.Bond(4, 1.0, true),
        Peridynamics.Bond(6, 1.0, true),
    ]
    @test b1.system.n_neighbors == [1, 2, 2, 2, 2]
    @test b1.system.bond_ids == [1:1, 2:3, 4:5, 6:7, 8:9]

    @test b1.storage.position == position[:, 1:6]
    @test b1.storage.displacement == zeros(3, 5)
    @test b1.storage.velocity == [1. 1. 1. 1. 1.; 0 0 0 0 0; 0 0 0 0 0]
    @test b1.storage.velocity_half == zeros(3, 5)
    @test b1.storage.acceleration == zeros(3, 5)
    @test b1.storage.b_int == zeros(3, 5)
    @test b1.storage.b_ext == zeros(3, 5)
    @test b1.storage.damage ≈ [0, 0, 0, 0, 0.5]
    @test b1.storage.bond_active == [1, 1, 1, 1, 1, 1, 1, 1, 0]
    @test b1.storage.n_active_bonds == [1, 2, 2, 2, 1]

    @test length(b1.paramsetup.parameters) == 2

    for i in 1:5
        pp = Peridynamics.get_params(b1, i)
        @test pp isa Peridynamics.StandardPointParameters
        @test pp.δ ≈ 1.01
        @test pp.rho ≈ 1.0
        @test pp.E ≈ 1.0
        @test pp.nu ≈ 0.25
        @test pp.G ≈ 0.4
        @test pp.K ≈ 2/3
        @test pp.λ ≈ 0.4
        @test pp.μ ≈ 0.4
        @test pp.Gc ≈ 1.0
        @test pp.εc ≈ 0.9083405243909494
        @test pp.bc ≈ 3.670674528926223
    end

    pp6 = Peridynamics.get_params(b1, 6)
    @test pp6 isa Peridynamics.StandardPointParameters
    @test pp6.δ ≈ 2.01
    @test pp6.rho ≈ 2.0
    @test pp6.E ≈ 2.0
    @test pp6.nu ≈ 0.25
    @test pp6.G ≈ 0.8
    @test pp6.K ≈ 4/3
    @test pp6.λ ≈ 0.8
    @test pp6.μ ≈ 0.8
    @test pp6.Gc ≈ 2.0
    @test pp6.εc ≈ 0.6438895077385466
    @test pp6.bc ≈ 0.4680337155970272

    condhandler1 = b1.condhandler
    @test condhandler1.loc_point_sets[:a] == [1, 2, 3, 4, 5]
    @test condhandler1.loc_point_sets[:b] == []
    @test length(condhandler1.single_dim_bcs) == 2
    @test condhandler1.single_dim_bcs[1].fun(0) == 0
    @test condhandler1.single_dim_bcs[1].fun(1) == 1
    @test condhandler1.single_dim_bcs[1].field === :velocity_half
    @test condhandler1.single_dim_bcs[1].point_set === :a
    @test condhandler1.single_dim_bcs[1].dim == 0x01

    ch1 = b1.system.chunk_handler
    @test ch1.point_ids == [1, 2, 3, 4, 5, 6]
    @test ch1.loc_points == 1:5
    @test ch1.halo_points == [6]
    @test ch1.hidxs_by_src[2] == 6:6
    for i in 1:6
        @test ch1.localizer[i] == i
    end

    b2 = body_chunks[2]
    @test b2 isa Peridynamics.BodyChunk
    @test b2.mat == BBMaterial()
    @test b2.system isa Peridynamics.BondSystem
    b2_ids = [6, 7, 8, 9, 10, 4, 5]
    @test b2.system.position == position[:, b2_ids]
    @test b2.system.volume == volume[b2_ids]
    @test b2.system.bonds == [
        Peridynamics.Bond(6, 2.0, true),
        Peridynamics.Bond(7, 1.0, true),
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(3, 2.0, true),
        Peridynamics.Bond(7, 2.0, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, 2.0, true),
        Peridynamics.Bond(1, 2.0, true),
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(4, 1.0, true),
        Peridynamics.Bond(5, 2.0, true),
        Peridynamics.Bond(2, 2.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(5, 1.0, true),
        Peridynamics.Bond(3, 2.0, true),
        Peridynamics.Bond(4, 1.0, true),
    ]
    @test b2.system.n_neighbors == [4, 4, 4, 3, 2]
    @test b2.system.bond_ids == [1:4, 5:8, 9:12, 13:15, 16:17]

    @test b2.storage.position == position[:, b2_ids]
    @test b2.storage.displacement == zeros(3, 5)
    @test b2.storage.velocity == zeros(3, 5)
    @test b2.storage.velocity_half == zeros(3, 5)
    @test b2.storage.acceleration == zeros(3, 5)
    @test b2.storage.b_int == zeros(3, 5)
    @test b2.storage.b_ext == zeros(3, 5)
    @test b2.storage.damage ≈ [0.5, 0.25, 0, 0, 0]
    @test b2.storage.bond_active == [0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    @test b2.storage.n_active_bonds == [4-2, 4-1, 4, 3, 2]

    @test length(b2.paramsetup.parameters) == 2

    for i in 1:5
        pp = Peridynamics.get_params(b2, i)
        @test pp isa Peridynamics.StandardPointParameters
        @test pp.δ ≈ 2.01
        @test pp.rho ≈ 2.0
        @test pp.E ≈ 2.0
        @test pp.nu ≈ 0.25
        @test pp.G ≈ 0.8
        @test pp.K ≈ 4/3
        @test pp.λ ≈ 0.8
        @test pp.μ ≈ 0.8
        @test pp.Gc ≈ 2.0
        @test pp.εc ≈ 0.6438895077385466
        @test pp.bc ≈ 0.4680337155970272
    end

    for i in 6:7
        pp = Peridynamics.get_params(b2, i)
        @test pp isa Peridynamics.StandardPointParameters
        @test pp.δ ≈ 1.01
        @test pp.rho ≈ 1.0
        @test pp.E ≈ 1.0
        @test pp.nu ≈ 0.25
        @test pp.G ≈ 0.4
        @test pp.K ≈ 2/3
        @test pp.λ ≈ 0.4
        @test pp.μ ≈ 0.4
        @test pp.Gc ≈ 1.0
        @test pp.εc ≈ 0.9083405243909494
        @test pp.bc ≈ 3.670674528926223
    end

    condhandler2 = b2.condhandler
    @test condhandler2.loc_point_sets[:a] == []
    @test condhandler2.loc_point_sets[:b] == [1, 2, 3, 4, 5]
    @test length(condhandler2.single_dim_bcs) == 2
    @test condhandler2.single_dim_bcs[1].fun(0) == 0
    @test condhandler2.single_dim_bcs[1].fun(1) == 1
    @test condhandler2.single_dim_bcs[1].field === :velocity_half
    @test condhandler2.single_dim_bcs[1].point_set === :a
    @test condhandler2.single_dim_bcs[1].dim == 0x01

    ch2 = b2.system.chunk_handler
    @test ch2.point_ids == b2_ids
    @test ch2.loc_points == 6:10
    @test ch2.halo_points == [4, 5]
    @test ch2.hidxs_by_src[1] == 6:7
    @test ch2.localizer[6] == 1
    @test ch2.localizer[7] == 2
    @test ch2.localizer[8] == 3
    @test ch2.localizer[9] == 4
    @test ch2.localizer[10] == 5
    @test ch2.localizer[4] == 6
    @test ch2.localizer[5] == 7
end
