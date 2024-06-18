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
