using Peridynamics, Test

let
    vv = VelocityVerlet(steps=1000)
    @test vv.end_time == -1.0
    @test vv.n_steps == 1000
    @test vv.Δt == -1.0
    @test vv.safety_factor == 0.7
end

let
    vv = VelocityVerlet(time=1)
    @test vv.end_time == 1.0
    @test vv.n_steps == -1
    @test vv.Δt == -1.0
    @test vv.safety_factor == 0.7
end

let
    vv = (@test_logs (:warn,) VelocityVerlet(time=1, stepsize=0.1))
    @test vv.end_time == 1.0
    @test vv.n_steps == -1
    @test vv.Δt == 0.1
    @test vv.safety_factor == 0.7
end

let
    vv = (@test_logs (:warn,) VelocityVerlet(steps=10, stepsize=0.1))
    @test vv.end_time == -1.0
    @test vv.n_steps == 10
    @test vv.Δt == 0.1
    @test vv.safety_factor == 0.7
end

let
    @test_throws ArgumentError VelocityVerlet(time=1, steps=10)
end

let
    @test_throws ArgumentError VelocityVerlet()
end

let
    @test_throws ArgumentError VelocityVerlet(time=1.0, safety_factor=0)
end

let
    @test_throws ArgumentError VelocityVerlet(time=1.0, safety_factor=1)
end
