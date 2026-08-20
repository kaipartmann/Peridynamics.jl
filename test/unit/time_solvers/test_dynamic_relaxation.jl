@testitem "DynamicRelaxation wrong input" begin
    @test_throws ArgumentError DynamicRelaxation(steps=0)
    @test_throws ArgumentError DynamicRelaxation(steps=10, stepsize=0)
    @test_throws ArgumentError DynamicRelaxation(steps=10, damping_factor=0)
end

@testitem "dynamic_relaxation_check" begin
    dr = DynamicRelaxation(steps=1)
    dr.n_steps = -1
    msg = "`n_steps` of DynamicRelaxation smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.dynamic_relaxation_check(dr)

    dr = DynamicRelaxation(steps=1)
    dr.Δt = -1
    msg = "`Δt` of DynamicRelaxation smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.dynamic_relaxation_check(dr)

    dr = DynamicRelaxation(steps=1)
    dr.Λ = -1
    msg = "`Λ` of DynamicRelaxation smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.dynamic_relaxation_check(dr)
end

@testitem "_init_density_matrix" begin
    using Peridynamics
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 0.0 1.0
                0.0 0.0 1.0 0.0]
    volume = fill(1.0, 4)
    body = Body(BBMaterial(), position, volume)
    material!(body; horizon=1.5, rho=8e-6, E=210e3)
    point_set!(body, :a, 1:2)
    material!(body, :a; horizon=1.5, rho=7.8e-6, E=200e3)

    dr = DynamicRelaxation(steps=1)
    dh = Peridynamics.threads_data_handler(body, dr, 1)
    chunk = dh.chunks[1]
    Peridynamics._init_density_matrix!(chunk,dr,chunk.paramsetup)

    @test chunk.storage.density_matrix ≈  [4.0e6  4.0e6  4.2e6  4.2e6
                                           4.0e6  4.0e6  4.2e6  4.2e6
                                           4.0e6  4.0e6  4.2e6  4.2e6]
end

@testitem "show DynamicRelaxation" begin
    io = IOBuffer()

    dr = DynamicRelaxation(steps=1)

    show(IOContext(io, :compact=>true), MIME("text/plain"), dr)
    msg = String(take!(io))
    @test contains(msg, "DynamicRelaxation(n_steps=1, Δt=1.0, Λ=1.0)")

    show(IOContext(io, :compact=>false), MIME("text/plain"), dr)
    msg = String(take!(io))
    @test contains(msg, "DynamicRelaxation:\n  n_steps  1\n  Δt       1\n  Λ        1")
end

@testitem "init_density_matrix! and relaxation_timestep! for all data handlers" begin
    import Peridynamics: init_density_matrix!, relaxation_timestep!
    # Setup single body
    position = [0.0 1.0 0.0 0.0;
                0.0 0.0 0.0 1.0;
                0.0 0.0 1.0 0.0]
    volume = fill(1.0, 4)
    body = Body(BBMaterial(), position, volume)
    material!(body; horizon=1.5, rho=8e-6, E=210e3)
    point_set!(body, :a, 1:2)
    material!(body, :a; horizon=3, rho=7.8e-6, E=200e3)
    forcedensity_bc!((p,t) -> 1e5 * p[1], body, :a, :x)
    dr = DynamicRelaxation(steps=2)
    options = Peridynamics.JobOptions(body)

    # ThreadsBodyDataHandler
    dh_threads = Peridynamics.threads_data_handler(body, dr, 1)
    init_density_matrix!(dh_threads, dr)
    @test all(dh_threads.chunks[1].storage.density_matrix .>= 0)
    @test dh_threads.chunks[1].storage.position ≈ position
    relaxation_timestep!(dh_threads, options, dr.Δt, 1)
    relaxation_timestep!(dh_threads, options, dr.Δt, 2)
    # Check that the density matrix is updated correctly
    @test dh_threads.chunks[1].storage.position[1, 2] ≈ 1.05

    # ThreadsMultibodyDataHandler
    b2 = Body(BBMaterial(), position, volume)
    material!(b2; horizon=1.5, rho=8e-6, E=210e3)
    ms = MultibodySetup(:b1 => body, :b2 => b2)
    dh_multibody = Peridynamics.threads_data_handler(ms, dr, 1)
    options_multi = Peridynamics.JobOptions(ms)
    init_density_matrix!(dh_multibody, dr)
    for body_idx in Peridynamics.each_body_idx(dh_multibody)
        for chunk in Peridynamics.get_body_dh(dh_multibody, body_idx).chunks
            @test all(chunk.storage.density_matrix .>= 0)
            @test chunk.storage.position ≈ position
        end
    end
    relaxation_timestep!(dh_multibody, options_multi, dr.Δt, 1)
    relaxation_timestep!(dh_multibody, options_multi, dr.Δt, 2)
    # Check that the density matrix is updated correctly
    @test dh_multibody.body_dhs[1].chunks[1].storage.position[1, 2] ≈ 1.05

    # AbstractMPIBodyDataHandler
    dh_mpi = Peridynamics.mpi_data_handler(body, dr)
    init_density_matrix!(dh_mpi, dr)
    @test all(dh_mpi.chunk.storage.density_matrix .>= 0)
    @test dh_mpi.chunk.storage.position ≈ position
    relaxation_timestep!(dh_mpi, options, dr.Δt, 1)
    relaxation_timestep!(dh_mpi, options, dr.Δt, 2)
    # Check that the density matrix is updated correctly
    @test dh_mpi.chunk.storage.position[1, 2] ≈ 1.05
end

@testitem "DynamicRelaxation: required storage fields" begin
    fields = Peridynamics.req_point_data_fields_timesolver(DynamicRelaxation)
    @test :velocity_half_old in fields && :b_int_old in fields && :density_matrix in fields
    @test Peridynamics.req_bond_data_fields_timesolver(DynamicRelaxation) == ()
    @test Peridynamics.req_data_fields_timesolver(DynamicRelaxation) == ()
end

@testitem "calc_damping: the damping coefficient is clamped" setup=[Fixtures] begin
    # cn = 2 sqrt(cn1 / cn2) with cn1 = -Σ u² Δb_int / (Δt ρ v_half_old), cn2 = Σ u²; a state
    # with a huge negative force increment gives cn > 2, which is clamped to 1.9
    body = Fixtures.line10()
    velocity_bc!(Fixtures.f_one, body, :all_points, :x)
    solver = DynamicRelaxation(steps=1)
    dh = Fixtures.handler(body, solver)
    Peridynamics.init_density_matrix!(dh, solver) # solve! does this before the time loop
    chunk = dh.chunks[1]
    (; displacement, velocity_half, velocity_half_old, b_int, b_int_old) = chunk.storage
    @test all(>(0), chunk.storage.density_matrix)
    displacement .= 1.0
    velocity_half .= 1.0
    velocity_half_old .= 1.0
    b_int_old .= 0.0
    b_int .= -1e6
    @test Peridynamics.calc_damping(chunk, 1.0) == 1.9
    # a positive force increment damps with the Rayleigh quotient itself
    b_int .= -1e-6
    cn1 = sum(1e-6 / dm for dm in chunk.storage.density_matrix) # all 30 dofs contribute
    @test Peridynamics.calc_damping(chunk, 1.0) ≈ 2 * sqrt(cn1 / 30)
    # a zero displacement field has no damping
    displacement .= 0.0
    @test Peridynamics.calc_damping(chunk, 1.0) == 0.0
end

@testitem "calc_damping: one global coefficient, independent of the decomposition" setup=[Fixtures] begin
    # the damping is a Rayleigh quotient over all dofs of the body; the chunked reduction and
    # the MPI Allgather path must reproduce the single-chunk value bitwise
    function damped_handler(n_chunks)
        body = Fixtures.line10()
        velocity_bc!(Fixtures.f_one, body, :all_points, :x)
        solver = DynamicRelaxation(steps=1)
        dh = Fixtures.handler(body, solver; n_chunks)
        Peridynamics.init_density_matrix!(dh, solver)
        for chunk in dh.chunks
            (; position, displacement, velocity_half, velocity_half_old,
               b_int, b_int_old) = chunk.storage
            n_loc = Peridynamics.get_n_loc_points(chunk.system)
            for (i, gi) in enumerate(chunk.system.chunk_handler.loc_points)
                displacement[:, i] .= 1e-3 * gi
                velocity_half[:, i] .= 0.1 * gi
                velocity_half_old[:, i] .= 0.1 * gi
                b_int[:, i] .= -1e2 * gi
                b_int_old[:, i] .= -1e2 * gi + 10.0
            end
        end
        return dh
    end
    cn_serial = Peridynamics.calc_damping(damped_handler(1).chunks[1], 1.0)
    @test 0 < cn_serial < 2
    for n_chunks in (2, 3, 5)
        cn = Peridynamics.calc_damping(damped_handler(n_chunks).chunks, 1.0)
        @test cn === cn_serial
    end
end
