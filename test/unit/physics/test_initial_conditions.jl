# Initial conditions (`src/physics/initial_conditions.jl`): declaration, conflicts and the
# values they set on a chunk.

@testitem "velocity_ic!: constant and position dependent values" setup=[Fixtures, ConditionCase] begin
    body = ConditionCase.two_points()
    velocity_ic!(body, :a, :x, 1.234)
    velocity_ic!(p -> p[1] * 3.456, body, :a, :y)
    velocity_ic!(body, :b, 3, -1.0) # the dimension as an integer
    chunk = Fixtures.handler(body, VelocityVerlet(steps=1); init=false).chunks[1]
    @test chunk.storage.velocity ≈ [1.234 1.234; 0.0 3.456; 0.0 -1.0]
end

@testitem "velocity_ic!: conflicts and argument checks" setup=[Fixtures, ConditionCase] begin
    F = Fixtures
    body = ConditionCase.two_points()
    velocity_ic!(body, :a, :x, 1.234)
    @test_throws ArgumentError velocity_ic!(body, :a, :x, 1.234)
    @test_throws ArgumentError velocity_ic!(F.f_p, body, :a, :x)
    @test_throws ArgumentError velocity_ic!(F.f_p, body, :b, :x)
    velocity_ic!(F.f_p, body, :a, :y)
    @test_throws ArgumentError velocity_ic!(F.f_p, body, :a, :y)
    @test_throws ArgumentError velocity_ic!(body, :a, :y, 1.234)
    @test_throws ArgumentError velocity_ic!(body, :b, :y, 1.234)
    # the function must be f(p)
    @test_throws ArgumentError velocity_ic!(F.f_bad_k, body, :a, :z)
    @test_throws ArgumentError velocity_ic!(F.f_bad_ab, body, :a, :z)
    # unknown dimensions and sets
    @test_throws ArgumentError velocity_ic!(body, :a, :k, 1.0)
    @test_throws ArgumentError velocity_ic!(body, :a, 4, 1.0)
    @test_throws ErrorException velocity_ic!(body, :c, :z, 1.0)
end
