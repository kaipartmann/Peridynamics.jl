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
    bonds = [Peridynamics.Bond(2, 1.0, true),
             Peridynamics.Bond(3, 1.0, true),
             Peridynamics.Bond(4, 1.0, true),
             Peridynamics.Bond(1, 1.0, true),
             Peridynamics.Bond(3, √2, true),
             Peridynamics.Bond(4, √2, true)]

    ch = Peridynamics.get_chunk_handler(bonds, pd, 1)
    @test ch.point_ids == [1, 2, 3, 4]
    @test ch.loc_points == 1:2
    @test ch.halo_points == [3, 4]
    @test ch.hidxs_by_src[2] == 3:4
    @test ch.localizer[1] == 1
    @test ch.localizer[2] == 2
    @test ch.localizer[3] == 3
    @test ch.localizer[4] == 4

    ch = Peridynamics.get_chunk_handler(bonds, pd, 2)
    @test ch.point_ids == [3, 4, 2, 1]
    @test ch.loc_points == 3:4
    @test ch.halo_points == [2, 1]
    @test ch.hidxs_by_src[1] == 3:4
    @test ch.localizer[1] == 4
    @test ch.localizer[2] == 3
    @test ch.localizer[3] == 1
    @test ch.localizer[4] == 2

    pd = Peridynamics.PointDecomposition(Peridynamics.distribute_equally(4, 4))
    ch = Peridynamics.get_chunk_handler(bonds, pd, 1)
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

    ch = Peridynamics.get_chunk_handler(bonds, pd, 2)
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
    N = 4
    pd = Peridynamics.PointDecomposition(Peridynamics.distribute_equally(N, 2))
    bonds = [Peridynamics.Bond(2, 1.0, true),
             Peridynamics.Bond(3, 1.0, true),
             Peridynamics.Bond(4, 1.0, true),
             Peridynamics.Bond(1, 1.0, true),
             Peridynamics.Bond(3, √2, true),
             Peridynamics.Bond(4, √2, true)]
    ch1 = Peridynamics.get_chunk_handler(bonds, pd, 1)
    ch2 = Peridynamics.get_chunk_handler(bonds, pd, 2)
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
