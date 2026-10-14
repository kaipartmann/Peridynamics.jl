@testitem "find_localizer" setup=[Fixtures] begin
    point_ids_1 = [10, 20, 30, 40]
    localizer_1 = Peridynamics.find_localizer(point_ids_1)
    @test localizer_1[10] == 1
    @test localizer_1[20] == 2
    @test localizer_1[30] == 3
    @test localizer_1[40] == 4

    point_ids_2 = unique(rand(Fixtures.rng(), 1:200, 200))
    localizer_2 = Peridynamics.find_localizer(point_ids_2)
    for (li, gi) in enumerate(point_ids_2)
        @test localizer_2[gi] == li
    end
end

@testitem "localize!" setup=[Fixtures] begin
    point_ids_1 = [10, 20, 30, 40]
    localizer_1 = Peridynamics.find_localizer(point_ids_1)
    point_set_1 = [40, 30, 20, 10]
    Peridynamics.localize!(point_set_1, localizer_1)
    @test point_set_1 == [4, 3, 2, 1]

    Peridynamics.localize!(point_ids_1, localizer_1)
    @test point_ids_1 == [1, 2, 3, 4]

    point_ids_2 = unique(rand(Fixtures.rng(), 1:200, 200))
    localizer_2 = Peridynamics.find_localizer(point_ids_2)
    Peridynamics.localize!(point_ids_2, localizer_2)
    for i in eachindex(point_ids_2)
        @test point_ids_2[i] == i
    end
end

@testitem "localize" setup=[Fixtures] begin
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

@testitem "localized_point_sets" setup=[Fixtures] begin
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

@testitem "localize!(Vector{Int}, ...)" setup=[Fixtures] begin
    # change two neighbors
    neighbor = [100, 101]
    localizer = Dict(100 => 1, 101 => 2)
    Peridynamics.localize!(neighbor, localizer)
    @test neighbor == [1, 2]

    # do not change any neighbor
    neighbor = [100, 101]
    localizer = Dict(100 => 100, 101 => 101)
    Peridynamics.localize!(neighbor, localizer)
    @test neighbor == [100, 101]

    # key not found error
    neighbor = [100]
    localizer = Dict(2 => 1)
    @test_throws KeyError(100) Peridynamics.localize!(neighbor, localizer)
end

@testitem "ChunkHandler" setup=[Fixtures] begin
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

@testitem "get_loc_view" setup=[Fixtures] begin
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

    neighbor1, bond_length1, fail_permit1, n_neighbors1, bond_ids1, ch1 = Peridynamics.get_bond_data(body, pd, 1)
    neighbor2, bond_length2, fail_permit2, n_neighbors2, bond_ids2, ch2 = Peridynamics.get_bond_data(body, pd, 2)
    rng = Fixtures.rng()
    v_float = rand(rng, N)
    m_float = rand(rng, 3, N)
    v_int = rand(rng, Int, N)
    m_int = rand(rng, Int, 3, N)

    @test Peridynamics.get_loc_view(v_int, ch1) == @view v_int[1:2]
    @test Peridynamics.get_loc_view(m_int, ch1) == @view m_int[:, 1:2]
    @test Peridynamics.get_loc_view(v_float, ch1) == @view v_float[1:2]
    @test Peridynamics.get_loc_view(m_float, ch1) == @view m_float[:, 1:2]
    @test Peridynamics.get_loc_view(v_int, ch2) == @view v_int[1:2]
    @test Peridynamics.get_loc_view(m_int, ch2) == @view m_int[:, 1:2]
    @test Peridynamics.get_loc_view(v_float, ch2) == @view v_float[1:2]
    @test Peridynamics.get_loc_view(m_float, ch2) == @view m_float[:, 1:2]
end

@testitem "sort_halo_by_src!" setup=[Fixtures] begin
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

@testitem "ChunkHandler, 10 points, 2 chunks" setup=[Fixtures] begin
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
    _neighbor1, _bond_length1, _fail_permit1, _n_neighbors1 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    neighbor1, bond_length1, fail_permit1, n_neighbors1, bond_ids1, ch1 = Peridynamics.get_bond_data(body, pd, chunk_id)
    @test _neighbor1 == neighbor1
    @test _bond_length1 == bond_length1
    @test _fail_permit1 == fail_permit1
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
    _neighbor2, _bond_length2, _fail_permit2, _n_neighbors2 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    neighbor2, bond_length2, fail_permit2, n_neighbors2, bond_ids2, ch2 = Peridynamics.get_bond_data(body, pd, chunk_id)
    @test ch2.n_loc_points == 5
    @test ch2.point_ids == [6, 7, 8, 9, 10, 5]
    @test ch2.loc_points == 6:10
    @test ch2.halo_points == [5]
    @test keys(ch2.hidxs_by_src) == Set([1])
    @test ch2.hidxs_by_src[1] == 6:6
    @test keys(ch2.localizer) == Set(5:10)
end

@testitem "ChunkHandler, 10 points, 5 chunks" setup=[Fixtures] begin
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
    _neighbor1, _bond_length1, _fail_permit1, _n_neighbors1 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    neighbor1, bond_length1, fail_permit1, n_neighbors1, bond_ids1, ch1 = Peridynamics.get_bond_data(body, pd, chunk_id)
    __neighbor1 = [2, 1, 3]
    __bond_length1 = [1.0, 1.0, 1.0]
    __fail_permit1 = [false, false, false]
    @test _neighbor1 == __neighbor1
    @test _bond_length1 == __bond_length1
    @test _fail_permit1 == __fail_permit1
    Peridynamics.localize!(__neighbor1, ch1.localizer)
    @test _neighbor1 == neighbor1 == __neighbor1
    @test bond_length1 == __bond_length1
    @test fail_permit1 == __fail_permit1
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
    _neighbor2, _bond_length2, _fail_permit2, _n_neighbors2 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    neighbor2, bond_length2, fail_permit2, n_neighbors2, bond_ids2, ch2 = Peridynamics.get_bond_data(body, pd, chunk_id)
    __neighbor2 = [2, 4, 3, 5]
    __bond_length2 = [1.0, 1.0, 1.0, 1.0]
    __fail_permit2 = [false, false, false, false]
    @test _neighbor2 == __neighbor2
    @test _bond_length2 == __bond_length2
    @test _fail_permit2 == __fail_permit2
    Peridynamics.localize!(__neighbor2, ch2.localizer)
    @test neighbor2 == __neighbor2
    @test bond_length2 == __bond_length2
    @test fail_permit2 == __fail_permit2
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
    _neighbor3, _bond_length3, _fail_permit3, _n_neighbors3 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    neighbor3, bond_length3, fail_permit3, n_neighbors3, bond_ids3, ch3 = Peridynamics.get_bond_data(body, pd, chunk_id)
    __neighbor3 = [4, 6, 5, 7]
    __bond_length3 = [1.0, 1.0, 1.0, 1.0]
    __fail_permit3 = [false, false, false, false]
    @test _neighbor3 == __neighbor3
    @test _bond_length3 == __bond_length3
    @test _fail_permit3 == __fail_permit3
    Peridynamics.localize!(__neighbor3, ch3.localizer)
    @test neighbor3 == __neighbor3
    @test bond_length3 == __bond_length3
    @test fail_permit3 == __fail_permit3
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
    _neighbor4, _bond_length4, _fail_permit4, _n_neighbors4 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    neighbor4, bond_length4, fail_permit4, n_neighbors4, bond_ids4, ch4 = Peridynamics.get_bond_data(body, pd, chunk_id)
    __neighbor4 = [6, 8, 7, 9]
    __bond_length4 = [1.0, 1.0, 1.0, 1.0]
    __fail_permit4 = [false, false, false, false]
    @test _neighbor4 == __neighbor4
    @test _bond_length4 == __bond_length4
    @test _fail_permit4 == __fail_permit4
    Peridynamics.localize!(__neighbor4, ch4.localizer)
    @test neighbor4 == __neighbor4
    @test bond_length4 == __bond_length4
    @test fail_permit4 == __fail_permit4
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
    _neighbor5, _bond_length5, _fail_permit5, _n_neighbors5 = Peridynamics.find_bonds(body, pd.decomp[chunk_id])
    neighbor5, bond_length5, fail_permit5, n_neighbors5, bond_ids5, ch5 = Peridynamics.get_bond_data(body, pd, chunk_id)
    __neighbor5 = [8, 10, 9]
    __bond_length5 = [1.0, 1.0, 1.0]
    __fail_permit5 = [false, false, false]
    @test _neighbor5 == __neighbor5
    @test _bond_length5 == __bond_length5
    @test _fail_permit5 == __fail_permit5
    Peridynamics.localize!(__neighbor5, ch5.localizer)
    @test neighbor5 == __neighbor5
    @test bond_length5 == __bond_length5
    @test fail_permit5 == __fail_permit5
    @test _n_neighbors5 == n_neighbors5 == [2, 1]
    @test ch5.n_loc_points == 2
    @test ch5.point_ids == [9, 10, 8]
    @test ch5.loc_points == 9:10
    @test ch5.halo_points == [8]
    @test keys(ch5.hidxs_by_src) == Set([4])
    @test ch5.hidxs_by_src[4] == 3:3
    @test keys(ch5.localizer) == Set(8:10)
end

@testitem "DeviceChunkHandler: the point counts survive, the bookkeeping stays home" setup=[Fixtures] begin
    # a minimal stand-in for the array type of another backend, see `test_storage_fields.jl`
    struct WrappedArray{T,N} <: AbstractArray{T,N}
        a::Array{T,N}
    end
    Base.size(x::WrappedArray) = size(x.a)
    Base.getindex(x::WrappedArray, i...) = getindex(x.a, i...)
    Base.setindex!(x::WrappedArray, v, i...) = setindex!(x.a, v, i...)
    struct WrappedBackend end
    function Peridynamics.Adapt.adapt_storage(::WrappedBackend, a::Array{T,N}) where {T,N}
        return WrappedArray{T,N}(a)
    end

    body = Fixtures.tetra4()
    chunk = Fixtures.chunk(body; n_chunks=2, chunk_id=1)
    ch = chunk.system.chunk_handler
    @test ch isa Peridynamics.ChunkHandler

    dch = Peridynamics.Adapt.adapt(WrappedBackend(), ch)
    @test dch isa Peridynamics.DeviceChunkHandler
    @test Peridynamics.get_n_loc_points(dch) == Peridynamics.get_n_loc_points(ch)
    @test Peridynamics.get_n_points(dch) == Peridynamics.get_n_points(ch)
    @test Peridynamics.each_point_idx(dch) == Peridynamics.each_point_idx(ch)

    # the local view of a field works on both, and it is the local points of the chunk
    a = collect(1.0:Peridynamics.get_n_points(ch))
    @test Peridynamics.get_loc_view(a, dch) == Peridynamics.get_loc_view(a, ch)

    # a target that moves no array hands back the chunk handler it was given, so that the
    # point ids, the halo bookkeeping and the localizer are not lost by an empty adapt
    @test Peridynamics.Adapt.adapt(Array, ch) === ch
end
