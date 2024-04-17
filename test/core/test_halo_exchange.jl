@testitem "find_halo_exchanges BBMaterial" begin
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

    loc_to_halo_exs, halo_to_loc_exs = Peridynamics.find_halo_exchanges(body_chunks)
    @test loc_to_halo_exs[1][1].src_chunk_id == 2
    @test loc_to_halo_exs[1][1].dest_chunk_id == 1
    @test loc_to_halo_exs[1][1].src_idxs == [1, 2]
    @test loc_to_halo_exs[1][1].dest_idxs == [3, 4]
    @test loc_to_halo_exs[2][1].src_chunk_id == 1
    @test loc_to_halo_exs[2][1].dest_chunk_id == 2
    @test loc_to_halo_exs[2][1].src_idxs == [1, 2]
    @test loc_to_halo_exs[2][1].dest_idxs == [3, 4]
    @test isempty(halo_to_loc_exs[1])
    @test isempty(halo_to_loc_exs[2])

    point_decomp = Peridynamics.PointDecomposition(body, 4)
    body_chunks = Peridynamics.chop_body_threads(body, ts, point_decomp, Val{1}())

    loc_to_halo_exs, halo_to_loc_exs = Peridynamics.find_halo_exchanges(body_chunks)

    he1 = sort(loc_to_halo_exs[1]; by=x -> x.src_chunk_id)
    @test he1[1].src_chunk_id == 2
    @test he1[1].dest_chunk_id == 1
    @test he1[1].src_idxs == [1]
    @test he1[1].dest_idxs == [2]
    @test he1[2].src_chunk_id == 3
    @test he1[2].dest_chunk_id == 1
    @test he1[2].src_idxs == [1]
    @test he1[2].dest_idxs == [3]
    @test he1[3].src_chunk_id == 4
    @test he1[3].dest_chunk_id == 1
    @test he1[3].src_idxs == [1]
    @test he1[3].dest_idxs == [4]

    he2 = sort(loc_to_halo_exs[2]; by=x -> x.src_chunk_id)
    @test he2[1].src_chunk_id == 1
    @test he2[1].dest_chunk_id == 2
    @test he2[1].src_idxs == [1]
    @test he2[1].dest_idxs == [2]
    @test he2[2].src_chunk_id == 3
    @test he2[2].dest_chunk_id == 2
    @test he2[2].src_idxs == [1]
    @test he2[2].dest_idxs == [3]
    @test he2[3].src_chunk_id == 4
    @test he2[3].dest_chunk_id == 2
    @test he2[3].src_idxs == [1]
    @test he2[3].dest_idxs == [4]

    he3 = sort(loc_to_halo_exs[3]; by=x -> x.src_chunk_id)
    @test he3[1].src_chunk_id == 1
    @test he3[1].dest_chunk_id == 3
    @test he3[1].src_idxs == [1]
    @test he3[1].dest_idxs == [2]
    @test he3[2].src_chunk_id == 2
    @test he3[2].dest_chunk_id == 3
    @test he3[2].src_idxs == [1]
    @test he3[2].dest_idxs == [3]
    @test he3[3].src_chunk_id == 4
    @test he3[3].dest_chunk_id == 3
    @test he3[3].src_idxs == [1]
    @test he3[3].dest_idxs == [4]

    he4 = sort(loc_to_halo_exs[4]; by=x -> x.src_chunk_id)
    @test he4[1].src_chunk_id == 1
    @test he4[1].dest_chunk_id == 4
    @test he4[1].src_idxs == [1]
    @test he4[1].dest_idxs == [2]
    @test he4[2].src_chunk_id == 2
    @test he4[2].dest_chunk_id == 4
    @test he4[2].src_idxs == [1]
    @test he4[2].dest_idxs == [3]
    @test he4[3].src_chunk_id == 3
    @test he4[3].dest_chunk_id == 4
    @test he4[3].src_idxs == [1]
    @test he4[3].dest_idxs == [4]
end

@testitem "exchange_read_fields! BBMaterial" begin
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
    b2 = tdh.chunks[2]
    @test b2 isa Peridynamics.BodyChunk
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

    randpos = rand(3, 4)
    tdh.chunks[2].storage.position .= randpos
    Peridynamics.exchange_read_fields!(tdh, 1)

    @test b1 isa Peridynamics.BodyChunk
    @test b1.storage.position[:,1:2] ≈ [0.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.position[:,3:4] ≈ randpos[:,1:2]
    @test b1.storage.displacement == zeros(3, 2)
    @test b1.storage.velocity == [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b1.storage.velocity_half == zeros(3, 2)
    @test b1.storage.acceleration == zeros(3, 2)
    @test b1.storage.b_int == zeros(3, 2)
    @test b1.storage.b_ext == zeros(3, 2)
    @test b1.storage.damage ≈ [2/3, 2/3]
    @test b1.storage.bond_active == [1, 0, 0, 1, 0, 0]
    @test b1.storage.n_active_bonds == [1, 1]

    @test b2 isa Peridynamics.BodyChunk
    @test b2.storage.position ≈ randpos
    @test b2.storage.displacement == zeros(3, 2)
    @test b2.storage.velocity == [0.0 0.0; 0.0 0.0; 0.0 0.0]
    @test b2.storage.velocity_half == zeros(3, 2)
    @test b2.storage.acceleration == zeros(3, 2)
    @test b2.storage.b_int == zeros(3, 2)
    @test b2.storage.b_ext == zeros(3, 2)
    @test b2.storage.damage ≈ [2/3, 2/3]
    @test b2.storage.bond_active == [0, 0, 1, 0, 0, 1]
    @test b2.storage.n_active_bonds == [1, 1]
end
