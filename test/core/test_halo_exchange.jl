@testitem "find_read_exchanges!" begin
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
    velocity_bc!(t -> t, body, :a, :x)
    forcedensity_bc!(t -> t, body, :a, :x)
    precrack!(body, :a, :b)
    ts = VelocityVerlet(steps=10)

    point_decomp = Peridynamics.PointDecomposition(body, 2)
    body_chunks = Peridynamics.chop_body_threads(body, ts, point_decomp)

    halo_exchanges = Peridynamics.find_halo_exchanges(body_chunks)
    @test halo_exchanges[1].field === :position
    @test halo_exchanges[1].from_chunk_id == 2
    @test halo_exchanges[1].to_chunk_id == 1
    @test halo_exchanges[1].from_loc_idxs == [1, 2]
    @test halo_exchanges[1].to_loc_idxs == [3, 4]
    @test halo_exchanges[2].field === :position
    @test halo_exchanges[2].from_chunk_id == 1
    @test halo_exchanges[2].to_chunk_id == 2
    @test halo_exchanges[2].from_loc_idxs == [1, 2]
    @test halo_exchanges[2].to_loc_idxs == [3, 4]

    point_decomp = Peridynamics.PointDecomposition(body, 4)
    body_chunks = Peridynamics.chop_body_threads(body, ts, point_decomp)

    halo_exchanges = Peridynamics.find_halo_exchanges(body_chunks)

    idx1 = findall(x -> x.to_chunk_id == 1, halo_exchanges)
    he1 = sort(halo_exchanges[idx1]; by = x -> x.from_chunk_id)

    @test he1[1].field === :position
    @test he1[1].from_chunk_id == 2
    @test he1[1].to_chunk_id == 1
    @test he1[1].from_loc_idxs == [1]
    @test he1[1].to_loc_idxs == [2]

    @test he1[2].field === :position
    @test he1[2].from_chunk_id == 3
    @test he1[2].to_chunk_id == 1
    @test he1[2].from_loc_idxs == [1]
    @test he1[2].to_loc_idxs == [3]

    @test he1[3].field === :position
    @test he1[3].from_chunk_id == 4
    @test he1[3].to_chunk_id == 1
    @test he1[3].from_loc_idxs == [1]
    @test he1[3].to_loc_idxs == [4]

    idx2 = findall(x -> x.to_chunk_id == 2, halo_exchanges)
    he2 = sort(halo_exchanges[idx2]; by = x -> x.from_chunk_id)

    @test he2[1].field === :position
    @test he2[1].from_chunk_id == 1
    @test he2[1].to_chunk_id == 2
    @test he2[1].from_loc_idxs == [1]
    @test he2[1].to_loc_idxs == [2]

    @test he2[2].field === :position
    @test he2[2].from_chunk_id == 3
    @test he2[2].to_chunk_id == 2
    @test he2[2].from_loc_idxs == [1]
    @test he2[2].to_loc_idxs == [3]

    @test he2[3].field === :position
    @test he2[3].from_chunk_id == 4
    @test he2[3].to_chunk_id == 2
    @test he2[3].from_loc_idxs == [1]
    @test he2[3].to_loc_idxs == [4]


    idx3 = findall(x -> x.to_chunk_id == 3, halo_exchanges)
    he3 = sort(halo_exchanges[idx3]; by = x -> x.from_chunk_id)

    @test he3[1].field === :position
    @test he3[1].from_chunk_id == 1
    @test he3[1].to_chunk_id == 3
    @test he3[1].from_loc_idxs == [1]
    @test he3[1].to_loc_idxs == [2]

    @test he3[2].field === :position
    @test he3[2].from_chunk_id == 2
    @test he3[2].to_chunk_id == 3
    @test he3[2].from_loc_idxs == [1]
    @test he3[2].to_loc_idxs == [3]

    @test he3[3].field === :position
    @test he3[3].from_chunk_id == 4
    @test he3[3].to_chunk_id == 3
    @test he3[3].from_loc_idxs == [1]
    @test he3[3].to_loc_idxs == [4]

    idx4 = findall(x -> x.to_chunk_id == 4, halo_exchanges)
    he4 = sort(halo_exchanges[idx4]; by = x -> x.from_chunk_id)

    @test he4[1].field === :position
    @test he4[1].from_chunk_id == 1
    @test he4[1].to_chunk_id == 4
    @test he4[1].from_loc_idxs == [1]
    @test he4[1].to_loc_idxs == [2]

    @test he4[2].field === :position
    @test he4[2].from_chunk_id == 2
    @test he4[2].to_chunk_id == 4
    @test he4[2].from_loc_idxs == [1]
    @test he4[2].to_loc_idxs == [3]

    @test he4[3].field === :position
    @test he4[3].from_chunk_id == 3
    @test he4[3].to_chunk_id == 4
    @test he4[3].from_loc_idxs == [1]
    @test he4[3].to_loc_idxs == [4]

end
