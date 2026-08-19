@testitem "BodyChunk: a 4-point body in 2 chunks, enumerated" begin
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

end

@testitem "BodyChunk: a 10-point line in 2 chunks, enumerated" begin
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

@testitem "chop_body_threads: every decomposition of a body is consistent" setup=[Fixtures] begin
    # Whatever the number of chunks, the chunks together are the body: the local points
    # partition the body, the halo points are exactly the neighbors of local points in other
    # chunks, the bonds of a chunk are the bonds of its local points with the lengths of the
    # body, the conditions and parameters are localized, and the pre-crack damage of every
    # point equals the damage the single-chunk body computes.
    position = zeros(3, 10)
    position[1, :] .= 0:9
    volume = collect(1.1:0.1:2.0)
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=1.01, rho=1, E=1, Gc=1)
    point_set!(body, :a, 1:5)
    point_set!(body, :b, 6:10)
    material!(body, :b, horizon=2.01, rho=2, E=2, Gc=2)
    velocity_ic!(body, :a, :x, 1.0)
    velocity_bc!(Fixtures.f_t, body, :a, :x)
    forcedensity_bc!(Fixtures.f_t, body, :a, :x)
    precrack!(body, :a, :b)
    ts = VelocityVerlet(steps=10)

    # the single-chunk reference
    ps = Peridynamics.get_param_spec(body)
    reference = Peridynamics.chop_body_threads(body, ts, Peridynamics.PointDecomposition(body, 1), ps)[1]
    @test Peridynamics.get_n_loc_points(reference) == Peridynamics.get_n_points(reference) == 10
    ref_damage = reference.storage.damage
    # set `a` (horizon 1.01) sees one neighbor on each side, set `b` (horizon 2.01) two; the
    # pre-crack breaks the bonds 5-6, 6-4 and 7-5, so point 5 loses 1 of 2, point 6 2 of 4
    # and point 7 1 of 4 bonds
    @test ref_damage ≈ [0, 0, 0, 0, 0.5, 0.5, 0.25, 0, 0, 0]

    for n_chunks in (2, 3, 5, 10)
        @testset "$n_chunks chunks" begin
            pd = Peridynamics.PointDecomposition(body, n_chunks)
            chunks = Peridynamics.chop_body_threads(body, ts, pd, ps)
            @test length(chunks) == n_chunks
            loc_ids = Int[]
            for (chunk_id, chunk) in enumerate(chunks)
                @test chunk isa Peridynamics.BodyChunk
                (; system, storage, paramsetup) = chunk
                ch = system.chunk_handler
                n_loc = Peridynamics.get_n_loc_points(chunk)
                @test n_loc == length(pd.decomp[chunk_id])
                @test ch.loc_points == pd.decomp[chunk_id]
                @test ch.point_ids[1:n_loc] == ch.loc_points
                append!(loc_ids, ch.loc_points)
                # the halo: neighbors of local points that live in other chunks, nothing else
                halo = ch.point_ids[(n_loc + 1):end]
                @test halo == ch.halo_points
                @test isempty(intersect(halo, ch.loc_points))
                neighbors_of_loc = Set{Int}()
                for (li, i) in enumerate(ch.loc_points)
                    for bond_id in Peridynamics.each_bond_idx(system, li)
                        j = ch.point_ids[system.bonds[bond_id].neighbor]
                        push!(neighbors_of_loc, j)
                        # the bond length is the distance of the points of the body
                        @test system.bonds[bond_id].length ≈ abs(position[1, i] - position[1, j])
                    end
                end
                @test Set(halo) == setdiff(neighbors_of_loc, ch.loc_points)
                # the localizer maps the global ids of the chunk back to local ids
                @test all(ch.localizer[ch.point_ids[k]] == k for k in eachindex(ch.point_ids))
                # positions and volumes are those of the body, local points first
                @test system.position == position[:, ch.point_ids]
                @test system.volume == volume[ch.point_ids]
                # parameters are mapped per local point
                @test all(Peridynamics.get_params(chunk, li).δ == body.point_params[body.params_map[i]].δ
                          for (li, i) in enumerate(ch.loc_points))
                # the initial condition on set `a`
                @test storage.velocity[1, :] == [i in 1:5 ? 1.0 : 0.0 for i in ch.loc_points]
                @test iszero(storage.velocity[2:3, :])
                # the pre-crack damage of the local points equals the reference
                @test storage.damage ≈ ref_damage[ch.loc_points]
                @test storage.n_active_bonds == [count(b -> storage.bond_active[b],
                                                       Peridynamics.each_bond_idx(system, li))
                                                 for li in 1:n_loc]
            end
            @test sort(loc_ids) == 1:10
        end
    end
end
