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
