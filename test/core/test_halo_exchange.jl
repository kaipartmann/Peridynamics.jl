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

    param_spec = Peridynamics.get_param_spec(body)
    point_decomp = Peridynamics.PointDecomposition(body, 2)
    body_chunks = Peridynamics.chop_body_threads(body, ts, point_decomp, param_spec)

    lth_exs, htl_exs = Peridynamics.find_halo_exchanges(body_chunks)

    @test lth_exs[1][1].src_chunk_id == 2
    @test lth_exs[1][1].dest_chunk_id == 1
    @test lth_exs[1][1].src_idxs == [1, 2]
    @test lth_exs[1][1].dest_idxs == [3, 4]
    @test lth_exs[2][1].src_chunk_id == 1
    @test lth_exs[2][1].dest_chunk_id == 2
    @test lth_exs[2][1].src_idxs == [1, 2]
    @test lth_exs[2][1].dest_idxs == [3, 4]

    @test htl_exs[1][1].src_chunk_id == 2
    @test htl_exs[1][1].dest_chunk_id == 1
    @test htl_exs[1][1].src_idxs == [3, 4]
    @test htl_exs[1][1].dest_idxs == [1, 2]
    @test htl_exs[2][1].src_chunk_id == 1
    @test htl_exs[2][1].dest_chunk_id == 2
    @test htl_exs[2][1].src_idxs == [3, 4]
    @test htl_exs[2][1].dest_idxs == [1, 2]

    point_decomp = Peridynamics.PointDecomposition(body, 4)
    body_chunks = Peridynamics.chop_body_threads(body, ts, point_decomp, param_spec)

    lth_exs, htl_exs = Peridynamics.find_halo_exchanges(body_chunks)

    he1 = sort(lth_exs[1]; by=x -> x.src_chunk_id)
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

    he2 = sort(lth_exs[2]; by=x -> x.src_chunk_id)
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

    he3 = sort(lth_exs[3]; by=x -> x.src_chunk_id)
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

    he4 = sort(lth_exs[4]; by=x -> x.src_chunk_id)
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

@testitem "exchange_loc_to_halo! BBMaterial" begin
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
    dh = Peridynamics.threads_data_handler(body, ts, 2)

    b1 = dh.chunks[1]
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
    b2 = dh.chunks[2]
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
    dh.chunks[2].storage.position .= randpos
    Peridynamics.exchange_loc_to_halo!(dh, 1)

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

@testitem "exchange core functions" begin
    buf_vec = zeros(10)
    buf_mat = zeros(3, 10)
    src_vec = zeros(20)
    dest_vec = zeros(20)
    src_mat = zeros(3, 20)
    dest_mat = zeros(3, 20)
    src_idxs = [1:10;]
    dest_idxs = [11:20;]

    # initialize src with values
    src_vec[src_idxs] .= 1.0
    src_mat[:, src_idxs] .= 1.0

    # exchange vector
    @test iszero(dest_vec[dest_idxs])
    Peridynamics.exchange!(dest_vec, src_vec, dest_idxs, src_idxs)
    @test dest_vec[dest_idxs] ≈ src_vec[src_idxs]
    Peridynamics.exchange_add!(dest_vec, src_vec, dest_idxs, src_idxs)
    @test dest_vec[dest_idxs] ≈ 2 * src_vec[src_idxs]

    # exchange matrix
    @test iszero(dest_mat[:, dest_idxs])
    Peridynamics.exchange!(dest_mat, src_mat, dest_idxs, src_idxs)
    @test dest_mat[:, dest_idxs] ≈ src_mat[:, src_idxs]
    Peridynamics.exchange_add!(dest_mat, src_mat, dest_idxs, src_idxs)
    @test dest_mat[:, dest_idxs] ≈ 2 * src_mat[:, src_idxs]

    # reset values and test again with buffer
    dest_vec .= 0.0
    dest_mat .= 0.0

    # exchange vector with buffer
    @test iszero(buf_vec)
    Peridynamics.exchange_to_buf!(buf_vec, src_vec, src_idxs)
    @test buf_vec[src_idxs] ≈ src_vec[src_idxs]
    Peridynamics.exchange_from_buf!(dest_vec, buf_vec, dest_idxs)
    @test dest_vec[dest_idxs] ≈ src_vec[src_idxs]
    @test buf_vec ≈ src_vec[src_idxs]
    Peridynamics.exchange_from_buf_add!(dest_vec, buf_vec, dest_idxs)
    @test dest_vec[dest_idxs] ≈ 2 * src_vec[src_idxs]
    @test buf_vec ≈ src_vec[src_idxs]

    # exchange matrix with buffer
    @test iszero(buf_mat)
    Peridynamics.exchange_to_buf!(buf_mat, src_mat, src_idxs)
    @test buf_mat[:, src_idxs] ≈ src_mat[:, src_idxs]
    Peridynamics.exchange_from_buf!(dest_mat, buf_mat, dest_idxs)
    @test dest_mat[:, dest_idxs] ≈ src_mat[:, src_idxs]
    @test buf_mat ≈ src_mat[:, src_idxs]
    Peridynamics.exchange_from_buf_add!(dest_mat, buf_mat, dest_idxs)
    @test dest_mat[:, dest_idxs] ≈ 2 * src_mat[:, src_idxs]
    @test buf_mat ≈ src_mat[:, src_idxs]
end

@testitem "local-to-halo exchange multithreading" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = CMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, nu=0.25, Gc=1)
    ts = VelocityVerlet(steps=10)
    dh = Peridynamics.threads_data_handler(body, ts, 2)

    # test data
    randpos = rand(3, 4)
    randbint = rand(3, 4)
    dh.chunks[2].storage.position .= randpos
    dh.chunks[2].storage.b_int .= randbint

    # exchange local-to-halo
    Peridynamics.exchange_loc_to_halo!(dh, 1)
    Peridynamics.exchange_loc_to_halo!(chunk -> chunk.storage.b_int, dh, 1)

    @test dh.chunks[1].storage.position[:,3:4] ≈ randpos[:,1:2]
    @test dh.chunks[1].storage.b_int[:,3:4] ≈ randbint[:,1:2]
end

@testitem "halo-to-local exchange multithreading" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = CMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, nu=0.25, Gc=1)
    ts = VelocityVerlet(steps=10)
    dh = Peridynamics.threads_data_handler(body, ts, 2)

    # test data
    randbint = rand(3, 4)
    randpos = rand(3, 4)
    dh.chunks[2].storage.b_int .= randbint
    dh.chunks[2].storage.position .= randpos

    # check that dest chunk is not modified
    @test iszero(dh.chunks[1].storage.b_int)
    @test dh.chunks[1].storage.position ≈ position

    # exchange local-to-halo
    Peridynamics.exchange_halo_to_loc!(dh, 1)
    Peridynamics.exchange_halo_to_loc!(chunk -> chunk.storage.position, dh, 1)

    @test dh.chunks[1].storage.b_int[:,1:2] ≈ randbint[:,3:4]
    @test dh.chunks[1].storage.position[:,1:2] ≈ position[:,1:2] + randpos[:,3:4]
end

@testitem "local-to-halo exchange MPI" tags=[:mpi] begin
    mpi_cmd = """
    using Peridynamics, Test
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = CMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, nu=0.25, Gc=1)
    ts = VelocityVerlet(steps=10)
    dh = Peridynamics.mpi_data_handler(body, ts)

    # test data
    randpos = [0.1 0.2 0.3 0.4
               0.1 0.2 0.3 0.4
               0.1 0.2 0.3 0.4]
    randbint = [0.11 0.21 0.31 0.41
                0.12 0.22 0.32 0.42
                0.13 0.23 0.33 0.43]
    rank = Peridynamics.mpi_rank()
    if rank == 1
        dh.chunk.storage.position .= randpos
        dh.chunk.storage.b_int .= randbint
    end

    # exchange local-to-halo
    Peridynamics.exchange_loc_to_halo!(dh)
    Peridynamics.exchange_loc_to_halo!(chunk -> chunk.storage.b_int, dh)

    if rank == 0
        @test dh.chunk.storage.position[:,3:4] ≈ randpos[:,1:2]
        @test dh.chunk.storage.b_int[:,3:4] ≈ randbint[:,1:2]
    end
    """
    mpiexec = Peridynamics.MPI.mpiexec()
    jlcmd = Base.julia_cmd()
    pdir = normpath(joinpath(@__DIR__, "..", ".."))
    cmd = `$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd)`
    @test success(cmd) # does not print anything
    # for debugging use the run command:
    # run(cmd)
end

@testitem "halo-to-local exchange MPI" tags=[:mpi] begin
    mpi_cmd = """
    using Peridynamics, Test
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = CMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, nu=0.25, Gc=1)
    ts = VelocityVerlet(steps=10)
    dh = Peridynamics.mpi_data_handler(body, ts)

    # test data
    randpos = [0.1 0.2 0.3 0.4
               0.1 0.2 0.3 0.4
               0.1 0.2 0.3 0.4]
    randbint = [0.11 0.21 0.31 0.41
                0.12 0.22 0.32 0.42
                0.13 0.23 0.33 0.43]
    rank = Peridynamics.mpi_rank()
    if rank == 1
        dh.chunk.storage.position .= randpos
        dh.chunk.storage.b_int .= randbint
    end

    # check that dest chunk is not modified
    if rank == 0
        @test iszero(dh.chunk.storage.b_int)
        @test dh.chunk.storage.position ≈ position
    end

    # exchange local-to-halo
    Peridynamics.exchange_halo_to_loc!(dh)
    Peridynamics.exchange_halo_to_loc!(chunk -> chunk.storage.position, dh)

    if rank == 0
        @test dh.chunk.storage.b_int[:,1:2] ≈ randbint[:,3:4]
        @test dh.chunk.storage.position[:,1:2] ≈ position[:,1:2] + randpos[:,3:4]
    end
    """
    mpiexec = Peridynamics.MPI.mpiexec()
    jlcmd = Base.julia_cmd()
    pdir = normpath(joinpath(@__DIR__, "..", ".."))
    cmd = `$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd)`
    @test success(cmd) # does not print anything
    # for debugging use the run command:
    # run(cmd)
end
