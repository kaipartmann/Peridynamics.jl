@testitem "find_localizer" begin
    point_ids_1 = [10, 20, 30, 40]
    localizer_1 = Peridynamics.find_localizer(point_ids_1)
    @test localizer_1[10] == 1
    @test localizer_1[20] == 2
    @test localizer_1[30] == 3
    @test localizer_1[40] == 4

    point_ids_2 = unique(rand(1:200, 200))
    localizer_2 = Peridynamics.find_localizer(point_ids_2)
    for (li, gi) in enumerate(point_ids_2)
        @test localizer_2[gi] == li
    end
end

@testitem "localize!" begin
    point_ids_1 = [10, 20, 30, 40]
    localizer_1 = Peridynamics.find_localizer(point_ids_1)
    point_set_1 = [40, 30, 20, 10]
    Peridynamics.localize!(point_set_1, localizer_1)
    @test point_set_1 == [4, 3, 2, 1]

    Peridynamics.localize!(point_ids_1, localizer_1)
    @test point_ids_1 == [1, 2, 3, 4]

    point_ids_2 = unique(rand(1:200, 200))
    localizer_2 = Peridynamics.find_localizer(point_ids_2)
    Peridynamics.localize!(point_ids_2, localizer_2)
    for i in eachindex(point_ids_2)
        @test point_ids_2[i] == i
    end
end

@testitem "localize" begin
    point_ids = collect(101:200)
    loc_points = 101:200
    n_loc_points = length(loc_points)
    halo_points = Vector{Int}()
    hidxs_by_src = Dict{Int,UnitRange{Int}}()
    localizer = Peridynamics.find_localizer(point_ids)
    ch = Peridynamics.ChunkHandler(n_loc_points, point_ids, loc_points, halo_points,
                                   hidxs_by_src, localizer)

    point_set = [101, 110, 120, 210, 220]
    loc_point_set = Peridynamics.localize(point_set, ch)
    @test loc_point_set == [1, 10, 20]
end

@testitem "localized_point_sets" begin
    point_ids = collect(101:200)
    loc_points = 101:200
    n_loc_points = length(loc_points)
    halo_points = Vector{Int}()
    hidxs_by_src = Dict{Int,UnitRange{Int}}()
    localizer = Peridynamics.find_localizer(point_ids)
    ch = Peridynamics.ChunkHandler(n_loc_points, point_ids, loc_points, halo_points,
                                   hidxs_by_src, localizer)

    point_sets = Dict(:a => [101, 110, 120, 210, 220], :b => [1, 2, 3])
    loc_point_sets = Peridynamics.localized_point_sets(point_sets, ch)
    @test loc_point_sets[:a] == [1, 10, 20]
    @test loc_point_sets[:b] == Vector{Int}()
end

@testitem "localize!(Vector{Bonds}, ...)" begin
    # change two bonds
    bonds = [Peridynamics.Bond(100, 1.0, true), Peridynamics.Bond(101, 30.0, false)]
    localizer = Dict(100 => 1, 101 => 2)
    Peridynamics.localize!(bonds, localizer)
    @test bonds == [Peridynamics.Bond(1, 1.0, true), Peridynamics.Bond(2, 30.0, false)]

    # do not change any bond
    bonds = [Peridynamics.Bond(100, 1.0, true), Peridynamics.Bond(101, 30.0, false)]
    localizer = Dict(100 => 100, 101 => 101)
    Peridynamics.localize!(bonds, localizer)
    @test bonds == [Peridynamics.Bond(100, 1.0, true), Peridynamics.Bond(101, 30.0, false)]

    # key not found error
    bonds = [Peridynamics.Bond(100, 1.0, true)]
    localizer = Dict(2 => 1)
    @test_throws KeyError(100) Peridynamics.localize!(bonds, localizer)
end

@testitem "ChunkHandler" begin
    pd = Peridynamics.PointDecomposition(Peridynamics.distribute_equally(4, 2))
    # bonds = [Peridynamics.Bond(2, 1.0, true),
    #          Peridynamics.Bond(3, 1.0, true),
    #          Peridynamics.Bond(4, 1.0, true),
    #          Peridynamics.Bond(1, 1.0, true),
    #          Peridynamics.Bond(3, √2, true),
    #          Peridynamics.Bond(4, √2, true)]

    halo_points = [3, 4]
    ch = Peridynamics.ChunkHandler(pd, halo_points, 1)
    @test ch.point_ids == [1, 2, 3, 4]
    @test ch.loc_points == 1:2
    @test ch.halo_points == [3, 4]
    @test ch.hidxs_by_src[2] == 3:4
    @test ch.localizer[1] == 1
    @test ch.localizer[2] == 2
    @test ch.localizer[3] == 3
    @test ch.localizer[4] == 4

    halo_points = [2, 1]
    ch = Peridynamics.ChunkHandler(pd, halo_points, 2)
    @test ch.point_ids == [3, 4, 2, 1]
    @test ch.loc_points == 3:4
    @test ch.halo_points == [2, 1]
    @test ch.hidxs_by_src[1] == 3:4
    @test ch.localizer[1] == 4
    @test ch.localizer[2] == 3
    @test ch.localizer[3] == 1
    @test ch.localizer[4] == 2

    pd = Peridynamics.PointDecomposition(Peridynamics.distribute_equally(4, 4))

    halo_points = [2, 3, 4]
    ch = Peridynamics.ChunkHandler(pd, halo_points, 1)
    @test ch.point_ids == [1, 2, 3, 4]
    @test ch.loc_points == 1:1
    @test ch.halo_points == [2, 3, 4]
    @test ch.hidxs_by_src[2] == 2:2
    @test ch.hidxs_by_src[3] == 3:3
    @test ch.hidxs_by_src[4] == 4:4
    @test ch.localizer[1] == 1
    @test ch.localizer[2] == 2
    @test ch.localizer[3] == 3
    @test ch.localizer[4] == 4

    halo_points = [1, 3, 4]
    ch = Peridynamics.ChunkHandler(pd, halo_points, 2)
    @test ch.point_ids == [2, 1, 3, 4]
    @test ch.loc_points == 2:2
    @test ch.halo_points == [1, 3, 4]
    @test ch.hidxs_by_src[1] == 2:2
    @test ch.hidxs_by_src[3] == 3:3
    @test ch.hidxs_by_src[4] == 4:4
    @test ch.localizer[1] == 2
    @test ch.localizer[2] == 1
    @test ch.localizer[3] == 3
    @test ch.localizer[4] == 4
end

@testitem "get_loc_view" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    N = length(volume)
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2.01, rho=1, E=1)
    n_chunks = 2
    decomp = Peridynamics.distribute_equally(body.n_points, n_chunks)
    pd = Peridynamics.PointDecomposition(decomp)

    bonds1, n_neighbors1, bond_ids1, ch1 = Peridynamics.get_bond_data(body, pd, 1)
    bonds2, n_neighbors2, bond_ids2, ch2 = Peridynamics.get_bond_data(body, pd, 2)
    v_float = rand(N)
    m_float = rand(3, N)
    v_int = rand(Int, N)
    m_int = rand(Int, 3, N)

    @test Peridynamics.get_loc_view(v_int, ch1) == @view v_int[1:2]
    @test Peridynamics.get_loc_view(m_int, ch1) == @view m_int[:, 1:2]
    @test Peridynamics.get_loc_view(v_float, ch1) == @view v_float[1:2]
    @test Peridynamics.get_loc_view(m_float, ch1) == @view m_float[:, 1:2]
    @test Peridynamics.get_loc_view(v_int, ch2) == @view v_int[1:2]
    @test Peridynamics.get_loc_view(m_int, ch2) == @view m_int[:, 1:2]
    @test Peridynamics.get_loc_view(v_float, ch2) == @view v_float[1:2]
    @test Peridynamics.get_loc_view(m_float, ch2) == @view m_float[:, 1:2]
end

@testitem "sort_halo_by_src!" begin
    # 10 points, 2 chunks, all halo points from same chunk
    halo_points = [8, 7, 6] # unsorted
    point_src = Dict(
        1 => 1, 2 => 1, 3 => 1, 4 => 1, 5 => 1, # all points in chunk 1
        6 => 2, 7 => 2, 8 => 2, 9 => 2, 10 => 2, # all points in chunk 2
    )
    n_loc_points = 5
    hidxs_by_src = Peridynamics.sort_halo_by_src!(halo_points, point_src, n_loc_points)
    @test halo_points == [8, 7, 6] # all in same chunk, order does not change
    @test hidxs_by_src[2] == 6:8

    # 10 points, 5 chunks, halo points from different chunks
    halo_points = [8,7,5,3] # unsorted
    point_src = Dict(
        1 => 1, 2 => 1, # all points in chunk 1
        3 => 2, 4 => 2, # all points in chunk 2
        5 => 3, 6 => 3, # all points in chunk 3
        7 => 4, 8 => 4, # all points in chunk 4
        9 => 5, 10 => 5, # all points in chunk 5
    )
    n_loc_points = 2
    hidxs_by_src = Peridynamics.sort_halo_by_src!(halo_points, point_src, n_loc_points)
    @test halo_points == [3, 5, 8, 7] # sorted by chunk, not fully sorted by id
    @test keys(hidxs_by_src) == Set([2, 3, 4])
    @test hidxs_by_src[2] == 3:3
    @test hidxs_by_src[3] == 4:4
    @test hidxs_by_src[4] == 5:6
end

@testitem "ChunkHandler, 10 points, 2 chunks" begin
    position = [0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    volume = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=1.01, rho=1, E=1)
    n_chunks = 2
    decomp = Peridynamics.distribute_equally(body.n_points, n_chunks)
    pd = Peridynamics.PointDecomposition(decomp)

    # chunk 1
    chunk_id = 1
    _bonds1, _n_neighbors1 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    bonds1, n_neighbors1, bond_ids1, ch1 = Peridynamics.get_bond_data(body, pd, chunk_id)
    @test _bonds1 == bonds1
    @test _n_neighbors1 == n_neighbors1
    @test ch1.n_loc_points == 5
    @test ch1.point_ids == [1, 2, 3, 4, 5, 6]
    @test ch1.loc_points == 1:5
    @test ch1.halo_points == [6]
    @test keys(ch1.hidxs_by_src) == Set([2])
    @test ch1.hidxs_by_src[2] == 6:6
    @test keys(ch1.localizer) == Set(1:6)

    # chunk 2
    chunk_id = 2
    _bonds2, _n_neighbors2 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    bonds2, n_neighbors2, bond_ids2, ch2 = Peridynamics.get_bond_data(body, pd, chunk_id)
    @test ch2.n_loc_points == 5
    @test ch2.point_ids == [6, 7, 8, 9, 10, 5]
    @test ch2.loc_points == 6:10
    @test ch2.halo_points == [5]
    @test keys(ch2.hidxs_by_src) == Set([1])
    @test ch2.hidxs_by_src[1] == 6:6
    @test keys(ch2.localizer) == Set(5:10)
end

@testitem "ChunkHandler, 10 points, 5 chunks" begin
    position = [0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    volume = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=1.01, rho=1, E=1)
    n_chunks = 5
    decomp = Peridynamics.distribute_equally(body.n_points, n_chunks)
    pd = Peridynamics.PointDecomposition(decomp)

    # chunk 1
    chunk_id = 1
    _bonds1, _n_neighbors1 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    bonds1, n_neighbors1, bond_ids1, ch1 = Peridynamics.get_bond_data(body, pd, chunk_id)
    __bonds1 = [
        Peridynamics.Bond(2, 1.0, false)
        Peridynamics.Bond(1, 1.0, false)
        Peridynamics.Bond(3, 1.0, false)
    ]
    @test _bonds1 == __bonds1
    Peridynamics.localize!(__bonds1, ch1.localizer)
    @test _bonds1 == bonds1 == __bonds1
    @test _n_neighbors1 == n_neighbors1 == [1, 2]
    @test ch1.n_loc_points == 2
    @test ch1.point_ids == [1, 2, 3]
    @test ch1.loc_points == 1:2
    @test ch1.halo_points == [3]
    @test keys(ch1.hidxs_by_src) == Set([2])
    @test ch1.hidxs_by_src[2] == 3:3
    @test keys(ch1.localizer) == Set(1:3)

    # chunk 2
    chunk_id = 2
    _bonds2, _n_neighbors2 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    bonds2, n_neighbors2, bond_ids2, ch2 = Peridynamics.get_bond_data(body, pd, chunk_id)
    __bonds2 = [
        Peridynamics.Bond(2, 1.0, false)
        Peridynamics.Bond(4, 1.0, false)
        Peridynamics.Bond(3, 1.0, false)
        Peridynamics.Bond(5, 1.0, false)
    ]
    @test _bonds2 == __bonds2
    Peridynamics.localize!(__bonds2, ch2.localizer)
    @test bonds2 == __bonds2
    @test _n_neighbors2 == n_neighbors2 == [2, 2]
    @test ch2.n_loc_points == 2
    @test ch2.point_ids == [3, 4, 2, 5]
    @test ch2.loc_points == 3:4
    @test ch2.halo_points == [2, 5]
    @test keys(ch2.hidxs_by_src) == Set([1, 3])
    @test ch2.hidxs_by_src[1] == 3:3
    @test ch2.hidxs_by_src[3] == 4:4
    @test keys(ch2.localizer) == Set(2:5)

    # chunk 3
    chunk_id = 3
    _bonds3, _n_neighbors3 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    bonds3, n_neighbors3, bond_ids3, ch3 = Peridynamics.get_bond_data(body, pd, chunk_id)
    __bonds3 = [
        Peridynamics.Bond(4, 1.0, false)
        Peridynamics.Bond(6, 1.0, false)
        Peridynamics.Bond(5, 1.0, false)
        Peridynamics.Bond(7, 1.0, false)
    ]
    @test _bonds3 == __bonds3
    Peridynamics.localize!(__bonds3, ch3.localizer)
    @test bonds3 == __bonds3
    @test _n_neighbors3 == n_neighbors3 == [2, 2]
    @test ch3.n_loc_points == 2
    @test ch3.point_ids == [5, 6, 4, 7]
    @test ch3.loc_points == 5:6
    @test ch3.halo_points == [4, 7]
    @test keys(ch3.hidxs_by_src) == Set([2, 4])
    @test ch3.hidxs_by_src[2] == 3:3
    @test ch3.hidxs_by_src[4] == 4:4
    @test keys(ch3.localizer) == Set(4:7)

    # chunk 4
    chunk_id = 4
    _bonds4, _n_neighbors4 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    bonds4, n_neighbors4, bond_ids4, ch4 = Peridynamics.get_bond_data(body, pd, chunk_id)
    __bonds4 = [
        Peridynamics.Bond(6, 1.0, false)
        Peridynamics.Bond(8, 1.0, false)
        Peridynamics.Bond(7, 1.0, false)
        Peridynamics.Bond(9, 1.0, false)
    ]
    @test _bonds4 == __bonds4
    Peridynamics.localize!(__bonds4, ch4.localizer)
    @test bonds4 == __bonds4
    @test _n_neighbors4 == n_neighbors4 == [2, 2]
    @test ch4.n_loc_points == 2
    @test ch4.point_ids == [7, 8, 6, 9]
    @test ch4.loc_points == 7:8
    @test ch4.halo_points == [6, 9]
    @test keys(ch4.hidxs_by_src) == Set([3, 5])
    @test ch4.hidxs_by_src[3] == 3:3
    @test ch4.hidxs_by_src[5] == 4:4
    @test keys(ch4.localizer) == Set(6:9)

    # chunk 5
    chunk_id = 5
    _bonds5, _n_neighbors5 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    bonds5, n_neighbors5, bond_ids5, ch5 = Peridynamics.get_bond_data(body, pd, chunk_id)
    __bonds5 = [
        Peridynamics.Bond(8, 1.0, false)
        Peridynamics.Bond(10, 1.0, false)
        Peridynamics.Bond(9, 1.0, false)
    ]
    @test _bonds5 == __bonds5
    Peridynamics.localize!(__bonds5, ch5.localizer)
    @test bonds5 == __bonds5
    @test _n_neighbors5 == n_neighbors5 == [2, 1]
    @test ch5.n_loc_points == 2
    @test ch5.point_ids == [9, 10, 8]
    @test ch5.loc_points == 9:10
    @test ch5.halo_points == [8]
    @test keys(ch5.hidxs_by_src) == Set([4])
    @test ch5.hidxs_by_src[4] == 3:3
    @test keys(ch5.localizer) == Set(8:10)
end
