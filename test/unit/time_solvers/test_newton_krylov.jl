@testitem "NewtonKrylov wrong input" begin
    @test_throws ArgumentError NewtonKrylov(time=1, steps=10)
    @test_throws ArgumentError NewtonKrylov()
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=-0.1)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, maxiter=0)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, maxiter=-1)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, tol=0)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, tol=-1e-8)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, gmres_reltol=0)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, gmres_reltol=-1e-4)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, gmres_abstol=0)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, gmres_abstol=-1e-8)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, gmres_maxiter=0)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, gmres_maxiter=-10)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, gmres_restart=0)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, gmres_restart=-10)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, perturbation_scale=0)
    @test_throws ArgumentError NewtonKrylov(steps=10, stepsize=0.1, perturbation_scale=-1)
end

@testitem "NewtonKrylov steps" begin
    nr = NewtonKrylov(steps=100)
    @test nr.end_time == 100
    @test nr.n_steps == 100
    @test nr.Δt == 1.0
    @test nr.maxiter == 100
    @test nr.tol == 1e-4
    @test nr.perturbation_scale == 1.0
    @test nr.gmres_maxiter == 200
    @test nr.gmres_reltol == 1e-4
    @test nr.gmres_abstol == 1e-8
end

@testitem "NewtonKrylov time" begin
    nr = NewtonKrylov(time=1.0, stepsize=0.01)
    @test nr.end_time == 1.0
    @test nr.n_steps == 100
    @test nr.Δt == 0.01
    @test nr.maxiter == 100
    @test nr.tol == 1e-4
end

@testitem "NewtonKrylov custom parameters" begin
    nr = NewtonKrylov(steps=50, stepsize=0.02, maxiter=75, tol=1e-6,
                      perturbation_scale=2.0, gmres_maxiter=150, gmres_reltol=1e-5,
                      gmres_abstol=1e-9)
    @test nr.end_time == 1.0
    @test nr.n_steps == 50
    @test nr.Δt == 0.02
    @test nr.maxiter == 75
    @test nr.tol == 1e-6
    @test nr.perturbation_scale == 2.0
    @test nr.gmres_maxiter == 150
    @test nr.gmres_reltol == 1e-5
    @test nr.gmres_abstol == 1e-9
end

@testitem "newton_krylov_check" begin
    nr = NewtonKrylov(steps=10, stepsize=0.1)
    nr.end_time = -1
    msg = "`end_time` of NewtonKrylov smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.newton_krylov_check(nr)

    nr = NewtonKrylov(steps=10, stepsize=0.1)
    nr.n_steps = -1
    msg = "`n_steps` of NewtonKrylov smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.newton_krylov_check(nr)

    nr = NewtonKrylov(steps=10, stepsize=0.1)
    nr.Δt = -1
    msg = "`Δt` of NewtonKrylov smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.newton_krylov_check(nr)

    nr = NewtonKrylov(steps=10, stepsize=0.1)
    nr.maxiter = -1
    msg = "`maxiter` of NewtonKrylov smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.newton_krylov_check(nr)

    nr = NewtonKrylov(steps=10, stepsize=0.1)
    nr.tol = -1
    msg = "`tol` of NewtonKrylov smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.newton_krylov_check(nr)

    nr = NewtonKrylov(steps=10, stepsize=0.1)
    nr.perturbation_scale = 0
    msg = "`perturbation_scale` of NewtonKrylov must be larger than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.newton_krylov_check(nr)
end

@testitem "show NewtonKrylov" begin
    io = IOBuffer()
    nr = NewtonKrylov(steps=10, stepsize=0.1)
    show(IOContext(io, :compact=>true), MIME("text/plain"), nr)
    msg = String(take!(io))
    @test contains(msg, "NewtonKrylov(end_time=1.0, n_steps=10, Δt=0.1, maxiter=100")

    show(IOContext(io, :compact=>false), MIME("text/plain"), nr)
    msg = String(take!(io))
    @test contains(msg, "NewtonKrylov:")
    @test contains(msg, "end_time")
    @test contains(msg, "n_steps")
    @test contains(msg, "Δt")
    @test contains(msg, "maxiter")
end

@testitem "chop_body_threads NewtonKrylov forces single chunk" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)

    nr = NewtonKrylov(steps=10, stepsize=0.1)
    point_decomp = Peridynamics.PointDecomposition(body, 4)
    param_spec = Peridynamics.get_param_spec(body)

    # NewtonKrylov now supports multiple chunks (parallel JFNK)
    chunks = Peridynamics.chop_body_threads(body, nr, point_decomp, param_spec)
    @test length(chunks) == 4  # Should use all requested chunks
end

@testitem "init_time_solver! NewtonKrylov ThreadsBodyDataHandler" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)
    point_set!(body, :top, [2])
    displacement_bc!(p -> 0.1, body, :top, :y)

    nr = NewtonKrylov(steps=10, stepsize=0.1)

    n_chunks = 1
    dh = Peridynamics.threads_data_handler(body, nr, n_chunks)
    Peridynamics.init_time_solver!(nr, dh)

    # Check that the buffers have been initialized
    n_points = size(body.position, 2)
    n_dof = 3 * n_points
    @test length(nr.global_residual) == n_dof
    @test length(nr.global_Δu) == n_dof
    @test length(nr.mpi_dof_counts) == 0
    @test length(nr.mpi_displs) == 0
    @test length(nr.u_norm_sq_buf) == n_chunks
    @test length(nr.v_norm_sq_buf) == n_chunks

    # Also multiple chunks
    n_chunks = 2
    dh = Peridynamics.threads_data_handler(body, nr, n_chunks)
    Peridynamics.init_time_solver!(nr, dh)

    # Check that the buffers have been initialized
    n_points = size(body.position, 2)
    n_dof = 3 * n_points
    @test length(nr.global_residual) == n_dof
    @test length(nr.global_Δu) == n_dof
    @test length(nr.mpi_dof_counts) == 0
    @test length(nr.mpi_displs) == 0
    @test length(nr.u_norm_sq_buf) == n_chunks
    @test length(nr.v_norm_sq_buf) == n_chunks
end

@testitem "init_time_solver! NewtonKrylov boundary condition checks" begin
    # Test: Body with displacement BC should work
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)
    point_set!(body, :top, [2])
    forcedensity_bc!(t -> 0.1, body, :top, :y)

    nr = NewtonKrylov(steps=10, stepsize=0.1)
    dh = Peridynamics.threads_data_handler(body, nr, 1)

    # a force density condition that is not position dependent is rejected (inside the
    # threaded loop, hence the CompositeException)
    @test_throws CompositeException Peridynamics.init_time_solver!(nr, dh)
end

@testitem "init_time_solver! NewtonKrylov damage check" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1.0)  # With Gc, damage is allowed
    point_set!(body, :a, [1, 2])
    point_set!(body, :b, [3, 4])
    precrack!(body, :a, :b)  # Create precrack
    point_set!(body, :top, [2])
    displacement_bc!(p -> 0.1, body, :top, :y)

    nr = NewtonKrylov(steps=10, stepsize=0.1)
    dh = Peridynamics.threads_data_handler(body, nr, 1)

    # should error, since damage currently not supported
    @test_throws CompositeException Peridynamics.init_time_solver!(nr, dh)
end

@testitem "init_time_solver! NewtonKrylov MultibodySetup error" begin
    position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    body1 = Body(BBMaterial(), position, volume)
    material!(body1, horizon=2, rho=1, E=1)
    point_set!(body1, :top, [2])
    displacement_bc!(p -> 0.1, body1, :top, :y)

    body2 = Body(BBMaterial(), position, volume)
    material!(body2, horizon=2, rho=1, E=1)

    ms = MultibodySetup(:b1 => body1, :b2 => body2)
    nr = NewtonKrylov(steps=10, stepsize=0.1)
    dh = Peridynamics.threads_data_handler(ms, nr, 1)

    msg = "NewtonKrylov solver only implemented for single body setups!\n"
    @test_throws ArgumentError(msg) Peridynamics.init_time_solver!(nr, dh)
end

@testitem "init_field_solver: NewtonKrylov claims its fields with markers" begin
    import Peridynamics: FullField, EmptyField, HaloPoints, init_field_solver
    nr = NewtonKrylov(steps=10, stepsize=0.1)
    vv = VelocityVerlet(steps=10)
    struct NKMarkerSystem <: Peridynamics.AbstractSystem end
    system = NKMarkerSystem()

    # the internal force density is assembled over the halo as well
    for f in (:b_int, :b_ext, :b_int_copy)
        @test init_field_solver(nr, system, Val(f)) === FullField(HaloPoints())
    end
    # the working fields of the solver, empty for every other solver
    for f in (:displacement_copy, :residual, :temp_force, :Δu, :v_temp, :Jv_temp)
        @test init_field_solver(nr, system, Val(f)) === FullField()
        @test init_field_solver(vv, system, Val(f)) === EmptyField()
    end
    # the fields every solver needs are left to the storage declaration
    @test isnothing(init_field_solver(nr, system, Val(:position)))
    @test isnothing(init_field_solver(nr, system, Val(:displacement)))
end

@testitem "NewtonKrylov: the storage fields have the extent the solver needs" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)

    nr = NewtonKrylov(steps=10, stepsize=0.1)
    dh = Peridynamics.threads_data_handler(body, nr, 1)
    (; system, storage) = dh.chunks[1]
    n_loc = Peridynamics.get_n_loc_points(system)
    n_all = Peridynamics.get_n_points(system)
    n_dof = Peridynamics.get_n_loc_dof(system)

    @test storage.position == position
    @test storage.displacement == zeros(3, n_loc)
    @test storage.displacement_copy == zeros(3, n_loc)
    for f in (:b_int, :b_ext, :b_int_copy)
        @test getfield(storage, f) == zeros(3, n_all)
    end
    for f in (:residual, :temp_force, :Δu, :v_temp, :Jv_temp)
        @test getfield(storage, f) == zeros(n_dof)
    end

    # the Velocity-Verlet storage of the same body leaves the solver fields empty
    vv = VelocityVerlet(steps=10)
    storage_vv = Peridynamics.threads_data_handler(body, vv, 1).chunks[1].storage
    @test size(storage_vv.displacement_copy) == (0, 0)
    @test size(storage_vv.b_int_copy) == (0, 0)
    for f in (:residual, :temp_force, :Δu, :v_temp, :Jv_temp)
        @test size(getfield(storage_vv, f)) == (0,)
    end
end
@testitem "req_point_data_fields_timesolver NewtonKrylov" begin
    fields = Peridynamics.req_point_data_fields_timesolver(NewtonKrylov)
    @test :displacement_copy in fields
    @test :b_int_copy in fields
    @test :residual in fields
    @test :v_temp in fields  # JFNK buffer
    @test :Jv_temp in fields  # JFNK buffer
end

@testitem "req_bond_data_fields_timesolver NewtonKrylov" begin
    fields = Peridynamics.req_bond_data_fields_timesolver(NewtonKrylov)
    @test fields == ()
end

@testitem "req_data_fields_timesolver NewtonKrylov" begin
    fields = Peridynamics.req_data_fields_timesolver(NewtonKrylov)
    @test fields == ()
end

@testitem "update_position!" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)
    point_set!(body, :top, [2])
    displacement_bc!(p -> 0.1, body, :top, :y)

    nr = NewtonKrylov(steps=10, stepsize=0.1)
    dh = Peridynamics.threads_data_handler(body, nr, 1)
    Peridynamics.init_time_solver!(nr, dh)

    chunk = dh.chunks[1]
    storage = chunk.storage
    system = chunk.system
    constrained_dofs = chunk.condhandler.constrained_dofs

    # Set some displacement values
    storage.displacement[2, 2] = 0.05  # y-displacement of point 2

    # Update position
    Peridynamics.update_position!(storage, system, constrained_dofs)

    # Check that position was updated for constrained DOF
    @test storage.position[2, 2] ≈ system.position[2, 2] + 0.05
end

@testitem "calc_residual!" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)
    point_set!(body, :top, [2])
    displacement_bc!(p -> 0.1, body, :top, :y)

    nr = NewtonKrylov(steps=10, stepsize=0.1)
    dh = Peridynamics.threads_data_handler(body, nr, 1)
    Peridynamics.init_time_solver!(nr, dh)

    chunk = dh.chunks[1]

    # Set some force values
    chunk.storage.b_int[1, 1] = 100.0
    chunk.storage.b_ext[1, 1] = -50.0

    # Calculate residual
    Peridynamics.calc_residual!(chunk)

    # Check that residual was calculated correctly (force * volume)
    expected_residual = (100.0 - 50.0) * volume[1]
    @test chunk.storage.residual[1] ≈ expected_residual

    # Check that constrained DOFs have zero residual
    for dof in chunk.condhandler.constrained_dofs
        @test chunk.storage.residual[dof] == 0.0
    end
end

@testitem "get_residual_norm" begin
    using Peridynamics.LinearAlgebra

    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)
    point_set!(body, :top, [2])
    displacement_bc!(p -> 0.1, body, :top, :y) # some dummy BC

    nr = NewtonKrylov(steps=10, stepsize=0.1)
    dh = Peridynamics.threads_data_handler(body, nr, 1)
    Peridynamics.init_time_solver!(nr, dh)

    chunk = dh.chunks[1]

    # Set some residual values
    chunk.storage.residual .= [1.0:12.0;]

    # the norm over all chunks of the data handler
    @test Peridynamics.get_residual_norm(dh) ≈ norm(chunk.storage.residual)
    # with two chunks: the norm of the concatenated residuals
    dh2 = Peridynamics.threads_data_handler(body, nr, 2)
    Peridynamics.init_time_solver!(nr, dh2)
    dh2.chunks[1].storage.residual .= 3.0
    dh2.chunks[2].storage.residual .= 4.0
    n1 = length(dh2.chunks[1].storage.residual)
    n2 = length(dh2.chunks[2].storage.residual)
    @test Peridynamics.get_residual_norm(dh2) ≈ sqrt(9n1 + 16n2)
end

@testitem "NewtonKrylov throw maxiter" tags=[:simulation] begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.0, 1.0, 1.0, 1.0]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)
    point_set!(body, :top, [4])
    point_set!(body, :bottom, [1])
    displacement_bc!(p -> 0.1, body, :top, :y)
    displacement_bc!(p -> 0.0, body, :bottom, :x)
    displacement_bc!(p -> 0.0, body, :bottom, :y)
    displacement_bc!(p -> 0.0, body, :bottom, :z)
    nr = NewtonKrylov(steps=10, stepsize=0.1, maxiter=2, tol=1e-6)
    job = Job(body, nr)
    @test_throws ErrorException submit(job; quiet=true)
end

@testitem "NewtonKrylov: the rotated formulations are rejected" begin
    # the materials with a stress rotation history are meant for dynamic simulations; the
    # solver rejects them when it is initialized on the data handler
    for mat in (CRMaterial(), RKCRMaterial())
        position = [0.0 1.0 0.0 0.0
                    0.0 0.0 1.0 0.0
                    0.0 0.0 0.0 1.0]
        body = Body(mat, position, [1.0, 1.0, 1.0, 1.0])
        material!(body, horizon=2, rho=1, E=1, nu=0.25)
        point_set!(body, :top, [4])
        displacement_bc!(p -> 0.1, body, :top, :y)
        nr = NewtonKrylov(steps=1, maxiter=1)
        dh = Peridynamics.threads_data_handler(body, nr, 1)
        @test_throws CompositeException Peridynamics.init_time_solver!(nr, dh)
    end
end
