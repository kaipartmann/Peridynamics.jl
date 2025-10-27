@testitem "NewtonRaphson wrong input" begin
    using Peridynamics: NewtonRaphson
    @test_throws ArgumentError NewtonRaphson(time=1, steps=10)
    @test_throws ArgumentError NewtonRaphson()
    @test_throws ArgumentError NewtonRaphson(steps=10, stepsize=0)
    @test_throws ArgumentError NewtonRaphson(steps=10, stepsize=-0.1)
    @test_throws ArgumentError NewtonRaphson(steps=10, stepsize=0.1, maxiter=0)
    @test_throws ArgumentError NewtonRaphson(steps=10, stepsize=0.1, maxiter=-1)
    @test_throws ArgumentError NewtonRaphson(steps=10, stepsize=0.1, tol=0)
    @test_throws ArgumentError NewtonRaphson(steps=10, stepsize=0.1, tol=-1e-8)
    @test_throws ArgumentError NewtonRaphson(steps=10, stepsize=0.1, gmres_reltol=0)
    @test_throws ArgumentError NewtonRaphson(steps=10, stepsize=0.1, gmres_reltol=-1e-4)
    @test_throws ArgumentError NewtonRaphson(steps=10, stepsize=0.1, gmres_abstol=0)
    @test_throws ArgumentError NewtonRaphson(steps=10, stepsize=0.1, gmres_abstol=-1e-8)
end

@testitem "NewtonRaphson steps" begin
    using Peridynamics: NewtonRaphson
    nr = NewtonRaphson(steps=100)
    @test nr.end_time == 100
    @test nr.n_steps == 100
    @test nr.Δt == 1.0
    @test nr.maxiter == 100
    @test nr.tol == 1e-4
    @test nr.perturbation == -1
    @test nr.gmres_maxiter == -1
    @test nr.gmres_reltol == 1e-4
    @test nr.gmres_abstol == 1e-8
end

@testitem "NewtonRaphson time" begin
    using Peridynamics: NewtonRaphson
    nr = NewtonRaphson(time=1.0, stepsize=0.01)
    @test nr.end_time == 1.0
    @test nr.n_steps == 100
    @test nr.Δt == 0.01
    @test nr.maxiter == 100
    @test nr.tol == 1e-4
end

@testitem "NewtonRaphson custom parameters" begin
    using Peridynamics: NewtonRaphson
    nr = NewtonRaphson(steps=50, stepsize=0.02, maxiter=75, tol=1e-6,
                      perturbation=1e-8, gmres_maxiter=150, gmres_reltol=1e-5,
                      gmres_abstol=1e-9)
    @test nr.end_time == 1.0
    @test nr.n_steps == 50
    @test nr.Δt == 0.02
    @test nr.maxiter == 75
    @test nr.tol == 1e-6
    @test nr.perturbation == 1e-8
    @test nr.gmres_maxiter == 150
    @test nr.gmres_reltol == 1e-5
    @test nr.gmres_abstol == 1e-9
end

@testitem "newton_raphson_check" begin
    using Peridynamics: NewtonRaphson
    nr = NewtonRaphson(steps=10, stepsize=0.1)
    nr.end_time = -1
    msg = "`end_time` of NewtonRaphson smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.newton_raphson_check(nr)

    nr = NewtonRaphson(steps=10, stepsize=0.1)
    nr.n_steps = -1
    msg = "`n_steps` of NewtonRaphson smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.newton_raphson_check(nr)

    nr = NewtonRaphson(steps=10, stepsize=0.1)
    nr.Δt = -1
    msg = "`Δt` of NewtonRaphson smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.newton_raphson_check(nr)

    nr = NewtonRaphson(steps=10, stepsize=0.1)
    nr.maxiter = -1
    msg = "`maxiter` of NewtonRaphson smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.newton_raphson_check(nr)

    nr = NewtonRaphson(steps=10, stepsize=0.1)
    nr.tol = -1
    msg = "`tol` of NewtonRaphson smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.newton_raphson_check(nr)
end

@testitem "show NewtonRaphson" begin
    using Peridynamics: NewtonRaphson

    io = IOBuffer()

    nr = NewtonRaphson(steps=10, stepsize=0.1)

    show(IOContext(io, :compact=>true), MIME("text/plain"), nr)
    msg = String(take!(io))
    @test contains(msg, "NewtonRaphson(end_time=1.0, n_steps=10, Δt=0.1, maxiter=100")

    show(IOContext(io, :compact=>false), MIME("text/plain"), nr)
    msg = String(take!(io))
    @test contains(msg, "NewtonRaphson:")
    @test contains(msg, "end_time")
    @test contains(msg, "n_steps")
    @test contains(msg, "Δt")
    @test contains(msg, "maxiter")
end

@testitem "chop_body_threads NewtonRaphson forces single chunk" begin
    using Peridynamics: NewtonRaphson
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)

    nr = NewtonRaphson(steps=10, stepsize=0.1)
    point_decomp = Peridynamics.PointDecomposition(body, 4)
    param_spec = Peridynamics.get_param_spec(body)

    # Test that warning is issued when more than 1 chunk is requested
    @test_logs (:warn,) begin
        chunks = Peridynamics.chop_body_threads(body, nr, point_decomp, param_spec)
        @test length(chunks) == 1  # Should force single chunk
    end
end

@testitem "init_time_solver! NewtonRaphson ThreadsBodyDataHandler" begin
    using Peridynamics: NewtonRaphson, displacement_bc!

    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)
    point_set!(body, :top, [2])
    displacement_bc!(p -> 0.1, body, :top, :y)

    nr = NewtonRaphson(steps=10, stepsize=0.1)

    dh = Peridynamics.threads_data_handler(body, nr, 1)
    Peridynamics.init_time_solver!(nr, dh)

    # Check that perturbation was set (if it was negative initially)
    @test nr.perturbation > 0

    # Check that GMRES maxiter was set
    @test nr.gmres_maxiter > 0
    @test nr.gmres_maxiter == min(200, Peridynamics.get_n_dof(dh.chunks[1].system))
end

@testitem "init_time_solver! NewtonRaphson checks single chunk" begin
    using Peridynamics: NewtonRaphson, displacement_bc!

    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)
    point_set!(body, :top, [2])
    displacement_bc!(p -> 0.1, body, :top, :y)

    nr = NewtonRaphson(steps=10, stepsize=0.1)

    # Create data handler with forced multiple chunks (should fail in init)
    # We need to manually create a dh with more than 1 chunk to test the error
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=10), 2)

    @test_throws ArgumentError Peridynamics.init_time_solver!(nr, dh)
end

@testitem "init_time_solver! NewtonRaphson boundary condition checks" begin
    using Peridynamics: NewtonRaphson, displacement_bc!

    # Test 1: No position-dependent BCs (should fail)
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body1 = Body(BBMaterial(), position, volume)
    material!(body1, horizon=2, rho=1, E=1)

    nr = NewtonRaphson(steps=10, stepsize=0.1)
    dh1 = Peridynamics.threads_data_handler(body1, nr, 1)

    @test_throws ArgumentError Peridynamics.init_time_solver!(nr, dh1)

    # Test 2: Time-dependent BCs (should fail)
    body2 = Body(BBMaterial(), position, volume)
    material!(body2, horizon=2, rho=1, E=1)
    point_set!(body2, :top, [2])
    velocity_bc!(t -> 0.1, body2, :top, :y)  # time-dependent BC

    nr2 = NewtonRaphson(steps=10, stepsize=0.1)
    dh2 = Peridynamics.threads_data_handler(body2, nr2, 1)

    @test_throws ArgumentError Peridynamics.init_time_solver!(nr2, dh2)
end

@testitem "init_time_solver! NewtonRaphson damage check" begin
    using Peridynamics: NewtonRaphson, displacement_bc!

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

    nr = NewtonRaphson(steps=10, stepsize=0.1)
    dh = Peridynamics.threads_data_handler(body, nr, 1)

    @test_throws ArgumentError Peridynamics.init_time_solver!(nr, dh)
end

@testitem "init_time_solver! NewtonRaphson MultibodySetup error" begin
    using Peridynamics: NewtonRaphson, displacement_bc!

    position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    body1 = Body(BBMaterial(), position, volume)
    material!(body1, horizon=2, rho=1, E=1)
    point_set!(body1, :top, [2])
    displacement_bc!(p -> 0.1, body1, :top, :y)

    body2 = Body(BBMaterial(), position, volume)
    material!(body2, horizon=2, rho=1, E=1)

    ms = MultibodySetup(:b1 => body1, :b2 => body2)
    nr = NewtonRaphson(steps=10, stepsize=0.1)
    dh = Peridynamics.threads_data_handler(ms, nr, 1)

    msg = "NewtonRaphson solver only implemented for single body multithreading!\n"
    @test_throws ArgumentError(msg) Peridynamics.init_time_solver!(nr, dh)
end

@testitem "init_field_solver NewtonRaphson" begin
    using Peridynamics: NewtonRaphson

    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)

    nr = NewtonRaphson(steps=10, stepsize=0.1)
    dh = Peridynamics.threads_data_handler(body, nr, 1)
    chunk = dh.chunks[1]
    system = chunk.system
    storage = chunk.storage

    # Test position field initialization
    pos_field = Peridynamics.init_field_solver(nr, system, Val(:position))
    @test pos_field == position

    # Test displacement field initialization
    disp_field = Peridynamics.init_field_solver(nr, system, Val(:displacement))
    @test disp_field == zeros(3, Peridynamics.get_n_loc_points(system))

    # Test b_int field initialization
    b_int_field = Peridynamics.init_field_solver(nr, system, Val(:b_int))
    @test b_int_field == zeros(3, Peridynamics.get_n_points(system))

    # Test b_ext field initialization
    b_ext_field = Peridynamics.init_field_solver(nr, system, Val(:b_ext))
    @test b_ext_field == zeros(3, Peridynamics.get_n_points(system))

    # Test displacement_copy field initialization
    disp_copy_field = Peridynamics.init_field_solver(nr, system, Val(:displacement_copy))
    @test disp_copy_field == zeros(3, Peridynamics.get_n_loc_points(system))

    # Test b_int_copy field initialization
    b_int_copy_field = Peridynamics.init_field_solver(nr, system, Val(:b_int_copy))
    @test b_int_copy_field == zeros(3, Peridynamics.get_n_points(system))

    # Test residual field initialization
    n_dof = Peridynamics.get_n_dof(system)
    residual_field = Peridynamics.init_field_solver(nr, system, Val(:residual))
    @test residual_field == zeros(n_dof)

    # Test jacobian field initialization
    jacobian_field = Peridynamics.init_field_solver(nr, system, Val(:jacobian))
    @test jacobian_field == zeros(n_dof, n_dof)

    # Test temp_force_a field initialization
    temp_a_field = Peridynamics.init_field_solver(nr, system, Val(:temp_force_a))
    @test temp_a_field == zeros(n_dof)

    # Test temp_force_b field initialization
    temp_b_field = Peridynamics.init_field_solver(nr, system, Val(:temp_force_b))
    @test temp_b_field == zeros(n_dof)

    # Test Δu field initialization
    du_field = Peridynamics.init_field_solver(nr, system, Val(:Δu))
    @test du_field == zeros(n_dof)

    # Test affected_points field initialization
    affected_field = Peridynamics.init_field_solver(nr, system, Val(:affected_points))
    @test length(affected_field) == n_dof ÷ 3
    @test all(x -> isa(x, Vector{Int}), affected_field)
end

@testitem "init_field_solver non-NewtonRaphson returns empty arrays" begin
    using Peridynamics: NewtonRaphson

    position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)
    vv = VelocityVerlet(steps=10)

    dh = Peridynamics.threads_data_handler(body, vv, 1)
    chunk = dh.chunks[1]
    system = chunk.system

    # Test that non-NR solver returns empty arrays for NR-specific fields
    disp_copy = Peridynamics.init_field_solver(vv, system, Val(:displacement_copy))
    @test size(disp_copy) == (0, 0)

    b_int_copy = Peridynamics.init_field_solver(vv, system, Val(:b_int_copy))
    @test size(b_int_copy) == (0, 0)

    residual = Peridynamics.init_field_solver(vv, system, Val(:residual))
    @test length(residual) == 0

    jacobian = Peridynamics.init_field_solver(vv, system, Val(:jacobian))
    @test size(jacobian) == (0, 0)

    temp_a = Peridynamics.init_field_solver(vv, system, Val(:temp_force_a))
    @test length(temp_a) == 0

    temp_b = Peridynamics.init_field_solver(vv, system, Val(:temp_force_b))
    @test length(temp_b) == 0

    du = Peridynamics.init_field_solver(vv, system, Val(:Δu))
    @test length(du) == 0

    affected = Peridynamics.init_field_solver(vv, system, Val(:affected_points))
    @test length(affected) == 0
end

@testitem "req_point_data_fields_timesolver NewtonRaphson" begin
    using Peridynamics: NewtonRaphson

    fields = Peridynamics.req_point_data_fields_timesolver(NewtonRaphson)

    @test :displacement_copy in fields
    @test :b_int_copy in fields
    @test :residual in fields
    @test :jacobian in fields
end

@testitem "req_bond_data_fields_timesolver NewtonRaphson" begin
    using Peridynamics: NewtonRaphson

    fields = Peridynamics.req_bond_data_fields_timesolver(NewtonRaphson)
    @test fields == ()
end

@testitem "req_data_fields_timesolver NewtonRaphson" begin
    using Peridynamics: NewtonRaphson

    fields = Peridynamics.req_data_fields_timesolver(NewtonRaphson)
    @test fields == ()
end

@testitem "update_position!" begin
    using Peridynamics: NewtonRaphson, displacement_bc!

    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)
    point_set!(body, :top, [2])
    displacement_bc!(p -> 0.1, body, :top, :y)

    nr = NewtonRaphson(steps=10, stepsize=0.1)
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
    using Peridynamics: NewtonRaphson, displacement_bc!

    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)
    point_set!(body, :top, [2])
    displacement_bc!(p -> 0.1, body, :top, :y)

    nr = NewtonRaphson(steps=10, stepsize=0.1)
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
    using Peridynamics: NewtonRaphson, displacement_bc!
    using Peridynamics.LinearAlgebra

    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)
    point_set!(body, :top, [2])
    displacement_bc!(p -> 0.1, body, :top, :y)

    nr = NewtonRaphson(steps=10, stepsize=0.1)
    dh = Peridynamics.threads_data_handler(body, nr, 1)
    Peridynamics.init_time_solver!(nr, dh)

    chunk = dh.chunks[1]

    # Set some residual values
    chunk.storage.residual .= [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]

    # Get residual norm from chunk
    r_chunk = Peridynamics.get_residual_norm(chunk)
    @test r_chunk ≈ norm(chunk.storage.residual)

    # Get residual norm from data handler
    r_dh = Peridynamics.get_residual_norm(dh)
    @test r_dh ≈ r_chunk
end

@testitem "minimum_volume ThreadsBodyDataHandler" begin
    using Peridynamics: NewtonRaphson, displacement_bc!

    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1)
    point_set!(body, :top, [2])
    displacement_bc!(p -> 0.1, body, :top, :y)

    nr = NewtonRaphson(steps=10, stepsize=0.1)
    dh = Peridynamics.threads_data_handler(body, nr, 1)

    min_vol = Peridynamics.minimum_volume(dh)
    @test min_vol == minimum(volume)
end
