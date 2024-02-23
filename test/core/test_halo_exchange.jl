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
    body_chunks = Peridynamics.chop_body_threads(body, ts, point_decomp, Val{1}())

    halo_exchanges = Peridynamics.find_halo_exchanges(body_chunks)
    @test halo_exchanges[1].field === :position
    @test halo_exchanges[1].src_chunk_id == 2
    @test halo_exchanges[1].dest_chunk_id == 1
    @test halo_exchanges[1].src_idxs == [1, 2]
    @test halo_exchanges[1].dest_idxs == [3, 4]
    @test halo_exchanges[2].field === :position
    @test halo_exchanges[2].src_chunk_id == 1
    @test halo_exchanges[2].dest_chunk_id == 2
    @test halo_exchanges[2].src_idxs == [1, 2]
    @test halo_exchanges[2].dest_idxs == [3, 4]

    point_decomp = Peridynamics.PointDecomposition(body, 4)
    body_chunks = Peridynamics.chop_body_threads(body, ts, point_decomp, Val{1}())

    halo_exchanges = Peridynamics.find_halo_exchanges(body_chunks)

    idx1 = findall(x -> x.dest_chunk_id == 1, halo_exchanges)
    he1 = sort(halo_exchanges[idx1]; by=x -> x.src_chunk_id)

    @test he1[1].field === :position
    @test he1[1].src_chunk_id == 2
    @test he1[1].dest_chunk_id == 1
    @test he1[1].src_idxs == [1]
    @test he1[1].dest_idxs == [2]

    @test he1[2].field === :position
    @test he1[2].src_chunk_id == 3
    @test he1[2].dest_chunk_id == 1
    @test he1[2].src_idxs == [1]
    @test he1[2].dest_idxs == [3]

    @test he1[3].field === :position
    @test he1[3].src_chunk_id == 4
    @test he1[3].dest_chunk_id == 1
    @test he1[3].src_idxs == [1]
    @test he1[3].dest_idxs == [4]

    idx2 = findall(x -> x.dest_chunk_id == 2, halo_exchanges)
    he2 = sort(halo_exchanges[idx2]; by=x -> x.src_chunk_id)

    @test he2[1].field === :position
    @test he2[1].src_chunk_id == 1
    @test he2[1].dest_chunk_id == 2
    @test he2[1].src_idxs == [1]
    @test he2[1].dest_idxs == [2]

    @test he2[2].field === :position
    @test he2[2].src_chunk_id == 3
    @test he2[2].dest_chunk_id == 2
    @test he2[2].src_idxs == [1]
    @test he2[2].dest_idxs == [3]

    @test he2[3].field === :position
    @test he2[3].src_chunk_id == 4
    @test he2[3].dest_chunk_id == 2
    @test he2[3].src_idxs == [1]
    @test he2[3].dest_idxs == [4]

    idx3 = findall(x -> x.dest_chunk_id == 3, halo_exchanges)
    he3 = sort(halo_exchanges[idx3]; by=x -> x.src_chunk_id)

    @test he3[1].field === :position
    @test he3[1].src_chunk_id == 1
    @test he3[1].dest_chunk_id == 3
    @test he3[1].src_idxs == [1]
    @test he3[1].dest_idxs == [2]

    @test he3[2].field === :position
    @test he3[2].src_chunk_id == 2
    @test he3[2].dest_chunk_id == 3
    @test he3[2].src_idxs == [1]
    @test he3[2].dest_idxs == [3]

    @test he3[3].field === :position
    @test he3[3].src_chunk_id == 4
    @test he3[3].dest_chunk_id == 3
    @test he3[3].src_idxs == [1]
    @test he3[3].dest_idxs == [4]

    idx4 = findall(x -> x.dest_chunk_id == 4, halo_exchanges)
    he4 = sort(halo_exchanges[idx4]; by=x -> x.src_chunk_id)

    @test he4[1].field === :position
    @test he4[1].src_chunk_id == 1
    @test he4[1].dest_chunk_id == 4
    @test he4[1].src_idxs == [1]
    @test he4[1].dest_idxs == [2]

    @test he4[2].field === :position
    @test he4[2].src_chunk_id == 2
    @test he4[2].dest_chunk_id == 4
    @test he4[2].src_idxs == [1]
    @test he4[2].dest_idxs == [3]

    @test he4[3].field === :position
    @test he4[3].src_chunk_id == 3
    @test he4[3].dest_chunk_id == 4
    @test he4[3].src_idxs == [1]
    @test he4[3].dest_idxs == [4]
end

@testitem "halo_exchange!" begin
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
    tdh = Peridynamics.ThreadsDataHandler(body, ts, point_decomp)

    b1 = tdh.chunks[1]
    @test b1 isa Peridynamics.BodyChunk
    @test b1.store.position == position
    @test b1.store.displacement == zeros(3, 2)
    @test b1.store.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.store.velocity_half == zeros(3, 2)
    @test b1.store.acceleration == zeros(3, 2)
    @test b1.store.b_int == zeros(3, 2)
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
    @test b2.store.b_int == zeros(3, 2)
    @test b2.store.b_ext == zeros(3, 2)
    @test b2.store.damage ≈ [2/3, 2/3]
    @test b2.store.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.store.n_active_bonds == [1, 1]

    randpos = rand(3, 4)
    tdh.chunks[2].store.position .= randpos
    Peridynamics.halo_exchange!(tdh, 1)

    @test b1 isa Peridynamics.BodyChunk
    @test b1.store.position[:,1:2] ≈ [0.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.store.position[:,3:4] ≈ randpos[:,1:2]
    @test b1.store.displacement == zeros(3, 2)
    @test b1.store.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.store.velocity_half == zeros(3, 2)
    @test b1.store.acceleration == zeros(3, 2)
    @test b1.store.b_int == zeros(3, 2)
    @test b1.store.b_ext == zeros(3, 2)
    @test b1.store.damage ≈ [2/3, 2/3]
    @test b1.store.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.store.n_active_bonds == [1, 1]

    @test b2 isa Peridynamics.BodyChunk
    @test b2.store.position ≈ randpos
    @test b2.store.displacement == zeros(3, 2)
    @test b2.store.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.store.velocity_half == zeros(3, 2)
    @test b2.store.acceleration == zeros(3, 2)
    @test b2.store.b_int == zeros(3, 2)
    @test b2.store.b_ext == zeros(3, 2)
    @test b2.store.damage ≈ [2/3, 2/3]
    @test b2.store.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.store.n_active_bonds == [1, 1]
end
