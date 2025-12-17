@testitem "VelocityVerlet steps" begin
    vv = VelocityVerlet(steps=1000)
    @test vv.end_time == -1.0
    @test vv.n_steps == 1000
    @test vv.Δt == -1.0
    @test vv.safety_factor == 0.7
end

@testitem "VelocityVerlet time" begin
    vv = VelocityVerlet(time=1)
    @test vv.end_time == 1.0
    @test vv.n_steps == -1
    @test vv.Δt == -1.0
    @test vv.safety_factor == 0.7
end

@testitem "VelocityVerlet time & stepsize" begin
    vv = (@test_logs (:warn,) VelocityVerlet(time=1, stepsize=0.1))
    @test vv.end_time == 1.0
    @test vv.n_steps == -1
    @test vv.Δt == 0.1
    @test vv.safety_factor == 0.7
end

@testitem "VelocityVerlet steps & stepsize" begin
    vv = (@test_logs (:warn,) VelocityVerlet(steps=10, stepsize=0.1))
    @test vv.end_time == -1.0
    @test vv.n_steps == 10
    @test vv.Δt == 0.1
    @test vv.safety_factor == 0.7
end

@testitem "VelocityVerlet wrong input" begin
    @test_throws ArgumentError VelocityVerlet(time=1, steps=10)
    @test_throws ArgumentError VelocityVerlet()
    @test_throws ArgumentError VelocityVerlet(time=1.0, safety_factor=0)
    @test_throws ArgumentError VelocityVerlet(time=1.0, safety_factor=1)
end

@testitem "init_time_solver! ThreadsBodyDataHandler" begin
    pos = [0 1; 0 0; 0 0]
    Δx, δ = 1, 1.5
    E, nu, rho = 1, 0.25, 1
    vol = fill(Δx^3, 2)
    body = Body(BBMaterial(), pos, vol)
    material!(body, horizon=δ, rho=rho, E=E, Gc=1)

    vv = VelocityVerlet(steps=10)

    dh = Peridynamics.threads_data_handler(body, vv, 1)
    Peridynamics.init_time_solver!(vv, dh)

    bc = 18 * E / (3 * (1 - 2 * nu)) / (π * δ^4)
    Δt = 0.7 * sqrt(2 * rho / bc)
    @test vv.Δt ≈ Δt atol=eps()
    @test vv.end_time ≈ 10Δt atol=eps()
    @test vv.n_steps == 10
    @test vv.safety_factor ≈ 0.7

    vv = VelocityVerlet(time=11)

    dh = Peridynamics.threads_data_handler(body, vv, 1)
    Peridynamics.init_time_solver!(vv, dh)

    bc = 18 * E / (3 * (1 - 2 * nu)) / (π * δ^4)
    Δt = 0.7 * sqrt(2 * rho / bc)
    @test vv.Δt ≈ Δt atol=eps()
    @test vv.n_steps == 10
    @test vv.end_time ≈ 11
    @test vv.safety_factor ≈ 0.7
end

@testitem "velocity_verlet_check" begin
    vv = VelocityVerlet(steps=1)
    msg = "`end_time` of VelocityVerlet smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.velocity_verlet_check(vv)

    vv = VelocityVerlet(time=1)
    msg = "`n_steps` of VelocityVerlet smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.velocity_verlet_check(vv)

    vv = VelocityVerlet(steps=1)
    vv.end_time = 1
    msg = "`Δt` of VelocityVerlet smaller than zero!\n"
    @test_throws ErrorException(msg) Peridynamics.velocity_verlet_check(vv)

    vv = VelocityVerlet(steps=1)
    vv.end_time = 1
    vv.Δt = 1
    vv.safety_factor = 0
    msg = "`safety_factor` of VelocityVerlet has invalid value!\n"
    @test_throws ErrorException(msg) Peridynamics.velocity_verlet_check(vv)
    vv.safety_factor = 1
    @test_throws ErrorException(msg) Peridynamics.velocity_verlet_check(vv)

end

@testitem "show VelocityVerlet" begin
    io = IOBuffer()

    vv = VelocityVerlet(steps=1)

    show(IOContext(io, :compact=>true), MIME("text/plain"), vv)
    msg = String(take!(io))
    @test contains(msg, "VelocityVerlet(n_steps=1, safety_factor=0.7)")

    show(IOContext(io, :compact=>false), MIME("text/plain"), vv)
    msg = String(take!(io))
    @test contains(msg, "VelocityVerlet:\n  n_steps        1\n  safety_factor  0.7")
end

@testitem "update_vel_half!" begin
    pos = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    vol = [1.0, 1.0]
    body = Body(BBMaterial(), pos, vol)
    material!(body, horizon=1.5, rho=1, E=1, Gc=1)
    vv = VelocityVerlet(steps=1)
    dh = Peridynamics.threads_data_handler(body, vv, 1)
    chunk = dh.chunks[1]

    # Set initial conditions
    chunk.storage.velocity .= [1.0 2.0; 3.0 4.0; 5.0 6.0]
    chunk.storage.acceleration .= [0.1 0.2; 0.3 0.4; 0.5 0.6]
    chunk.storage.velocity_half .= 0.0

    Δt½ = 0.5
    Peridynamics.update_vel_half!(chunk, Δt½)

    # Check: velocity_half = velocity + acceleration * Δt½
    expected = [1.0 2.0; 3.0 4.0; 5.0 6.0] .+ [0.1 0.2; 0.3 0.4; 0.5 0.6] .* Δt½
    @test chunk.storage.velocity_half ≈ expected
end

@testitem "update_disp_and_pos!" begin
    pos = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    vol = [1.0, 1.0]
    body = Body(BBMaterial(), pos, vol)
    material!(body, horizon=1.5, rho=1, E=1, Gc=1)
    vv = VelocityVerlet(steps=1)
    dh = Peridynamics.threads_data_handler(body, vv, 1)
    chunk = dh.chunks[1]

    # Set initial conditions
    initial_displacement = [0.1 0.2; 0.3 0.4; 0.5 0.6]
    initial_position = copy(chunk.storage.position)
    chunk.storage.displacement .= initial_displacement
    chunk.storage.velocity_half .= [1.0 2.0; 3.0 4.0; 5.0 6.0]

    Δt = 0.1
    Peridynamics.update_disp_and_pos!(chunk, Δt)

    # Check: displacement += velocity_half * Δt
    #        position += velocity_half * Δt
    du = [1.0 2.0; 3.0 4.0; 5.0 6.0] .* Δt
    @test chunk.storage.displacement ≈ initial_displacement .+ du
    @test chunk.storage.position ≈ initial_position .+ du
end

@testitem "_update_acc! and _update_vel!" begin
    # Test the inline helper functions directly
    acceleration = zeros(6)
    velocity = zeros(6)
    velocity_half = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    b_int = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
    b_ext = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    rho = 2.0
    Δt½ = 0.5

    for dof in 1:6
        Peridynamics._update_acc!(acceleration, b_int, b_ext, rho, dof)
        Peridynamics._update_vel!(velocity, velocity_half, acceleration, Δt½, dof)
    end

    # Check: acceleration = (b_int + b_ext) / rho
    expected_acc = (b_int .+ b_ext) ./ rho
    @test acceleration ≈ expected_acc

    # Check: velocity = velocity_half + acceleration * Δt½
    expected_vel = velocity_half .+ expected_acc .* Δt½
    @test velocity ≈ expected_vel
end

@testitem "update_acc_and_vel! with uniform parameters" begin
    # Test with uniform parameters (single AbstractPointParameters)
    pos = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    vol = [1.0, 1.0]
    body = Body(BBMaterial(), pos, vol)
    material!(body, horizon=1.5, rho=2.0, E=1, Gc=1)
    vv = VelocityVerlet(steps=1)
    dh = Peridynamics.threads_data_handler(body, vv, 1)
    chunk = dh.chunks[1]

    # Set initial conditions
    chunk.storage.velocity_half .= [1.0 2.0; 3.0 4.0; 5.0 6.0]
    chunk.storage.b_int .= [10.0 20.0; 30.0 40.0; 50.0 60.0]
    chunk.storage.b_ext .= [1.0 2.0; 3.0 4.0; 5.0 6.0]
    chunk.storage.acceleration .= 0.0
    chunk.storage.velocity .= 0.0

    Δt½ = 0.5
    Peridynamics.update_acc_and_vel!(chunk, Δt½)

    # Check: acceleration = (b_int + b_ext) / rho
    rho = 2.0
    expected_acc = (chunk.storage.b_int .+ chunk.storage.b_ext) ./ rho
    @test chunk.storage.acceleration ≈ expected_acc

    # Check: velocity = velocity_half + acceleration * Δt½
    expected_vel = [1.0 2.0; 3.0 4.0; 5.0 6.0] .+ expected_acc .* Δt½
    @test chunk.storage.velocity ≈ expected_vel
end

@testitem "update_acc_and_vel! with parameter handler" begin
    # Test with AbstractParameterHandler (multiple point parameters)
    pos = [0.0 1.0 2.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    vol = [1.0, 1.0, 1.0]
    body = Body(BBMaterial(), pos, vol)
    material!(body, horizon=1.5, rho=1.0, E=1, Gc=1)
    # Add a second parameter set with different density
    point_set!(body, :set2, [3])
    material!(body, :set2, horizon=1.5, rho=3.0, E=1, Gc=1)

    vv = VelocityVerlet(steps=1)
    dh = Peridynamics.threads_data_handler(body, vv, 1)
    chunk = dh.chunks[1]

    # Set initial conditions
    chunk.storage.velocity_half .= [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    chunk.storage.b_int .= [10.0 20.0 30.0; 40.0 50.0 60.0; 70.0 80.0 90.0]
    chunk.storage.b_ext .= [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    chunk.storage.acceleration .= 0.0
    chunk.storage.velocity .= 0.0

    Δt½ = 0.5
    Peridynamics.update_acc_and_vel!(chunk, Δt½)

    # Check that accelerations are computed correctly for each point
    # Points 1 and 2 should use rho=1.0, point 3 should use rho=3.0
    for i in 1:3
        params = Peridynamics.get_params(chunk, i)
        for dim in 1:3
            dof = (i - 1) * 3 + dim
            expected_acc = (chunk.storage.b_int[dof] + chunk.storage.b_ext[dof]) / params.rho
            @test chunk.storage.acceleration[dof] ≈ expected_acc

            expected_vel = chunk.storage.velocity_half[dof] + expected_acc * Δt½
            @test chunk.storage.velocity[dof] ≈ expected_vel
        end
    end
end
