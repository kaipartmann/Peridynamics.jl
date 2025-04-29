@testitem "All conditions" begin
    position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=1.5, rho=1, E=210e9, Gc=1.0)
    point_set!(body, :a, 1:2)
    velocity_bc!(t -> t, body, :a, :x)
    velocity_bc!((p,t) -> p[1] * t, body, :a, :y)
    forcedensity_bc!(t -> t, body, :a, :x)
    forcedensity_bc!((p,t) -> p[1] * t, body, :a, :y)
    velocity_ic!(body, :a, :x, 1.234)
    velocity_ic!(p -> p[1] * 3.456, body, :a, :y)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; velocity_half, velocity, b_ext) = chunk.storage

    @test velocity ≈ [1.234 1.234; 0.0 3.456; 0.0 0.0]

    Peridynamics.apply_boundary_conditions!(chunk, 2.345)

    @test velocity_half ≈ [2.345 2.345; 0.0 2.345; 0.0 0.0]
    @test b_ext ≈ [2.345 2.345; 0.0 2.345; 0.0 0.0]
end

@testitem "Condition errors" begin
    position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=1.5, rho=1, E=210e9, Gc=1.0)
    point_set!(body, :a, 1:2)
    point_set!(body, :b, [2])

    # conflicts with existing conditions:
    velocity_bc!(t -> t, body, :a, :x)
    @test_throws ArgumentError velocity_bc!(t -> 2 * t, body, :a, :x)
    @test_throws ArgumentError velocity_bc!((p,t) -> p[1] * t, body, :a, :x)
    @test_throws ArgumentError velocity_bc!((p,t) -> p[1] * t, body, :b, :x)

    velocity_bc!((p,t) -> p[1] * t, body, :a, :y)
    @test_throws ArgumentError velocity_bc!((p,t) -> p[2] * t, body, :a, :y)
    @test_throws ArgumentError velocity_bc!(t -> t, body, :a, :y)
    @test_throws ArgumentError velocity_bc!(t -> t, body, :b, :y)

    forcedensity_bc!(t -> t, body, :a, :x)
    @test_throws ArgumentError forcedensity_bc!(t -> 2 * t, body, :a, :x)
    @test_throws ArgumentError forcedensity_bc!((p,t) -> p[1] * t, body, :a, :x)
    @test_throws ArgumentError forcedensity_bc!((p,t) -> p[1] * t, body, :b, :x)

    forcedensity_bc!((p,t) -> p[1] * t, body, :a, :y)
    @test_throws ArgumentError forcedensity_bc!((p,t) -> p[2] * t, body, :a, :y)
    @test_throws ArgumentError forcedensity_bc!(t -> t, body, :a, :y)
    @test_throws ArgumentError forcedensity_bc!(t -> t, body, :b, :y)

    velocity_ic!(body, :a, :x, 1.234)
    @test_throws ArgumentError velocity_ic!(body, :a, :x, 1.234)
    @test_throws ArgumentError velocity_ic!(p -> 3.456, body, :a, :x)
    @test_throws ArgumentError velocity_ic!(p -> 3.456, body, :b, :x)

    velocity_ic!(p -> 3.456, body, :a, :y)
    @test_throws ArgumentError velocity_ic!(p -> 3.456, body, :a, :y)
    @test_throws ArgumentError velocity_ic!(body, :a, :y, 1.234)
    @test_throws ArgumentError velocity_ic!(body, :b, :y, 1.234)

    # Wrong condition function arguments
    @test_throws ArgumentError velocity_bc!((a, b) -> a * b, body, :a, :z)
    @test_throws ArgumentError velocity_bc!((k, t, u) -> k * t * u, body, :a, :z)
    @test_throws ArgumentError forcedensity_bc!((a, b) -> a * b, body, :a, :z)
    @test_throws ArgumentError forcedensity_bc!((k, t, u) -> k * t * u, body, :a, :z)
    @test_throws ArgumentError velocity_ic!(k -> 3.456, body, :a, :z)
    @test_throws ArgumentError velocity_ic!((a, b) -> 3.456, body, :a, :z)

    # unknown symbols
    @test_throws ArgumentError velocity_bc!(t -> 1, body, :a, :k)
    @test_throws ArgumentError velocity_bc!(t -> 1, body, :a, 4)
end

@testitem "Condition errors DataBCs" begin
    position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=1.5, rho=1, E=210e9, Gc=1.0)
    point_set!(body, :a, 1:2)
    point_set!(body, :b, [2])

    # data has wrong dimensions
    data = rand(2, 3) # not matching the number of points
    @test_throws ArgumentError Peridynamics.velocity_databc!(body, data, :a, [1, 2])
    data = rand(3, 2) # not matching the number of dimensions
    @test_throws ArgumentError Peridynamics.velocity_databc!(body, data, :a, [1, 2])

    # should work for just the x-dimension
    data = rand(1, 2)
    Peridynamics.velocity_databc!(body, data, :a, [1])

    # conflicts with existing conditions:
    @test_throws ArgumentError velocity_bc!(t -> 2 * t, body, :a, :x)

    data = rand(1, 2)
    Peridynamics.velocity_databc!(body, data, :a, [:y])

    @test_throws ArgumentError velocity_bc!(t -> t, body, :a, :y)
    @test_throws ArgumentError velocity_bc!(t -> t, body, :b, :y)

    data = rand(1, 2)
    Peridynamics.forcedensity_databc!(body, data, :a, [:x])
    @test_throws ArgumentError Peridynamics.forcedensity_databc!(body, 2 * data, :a, [:x])
    @test_throws ArgumentError forcedensity_bc!((p,t) -> p[1] * t, body, :a, :x)
    @test_throws ArgumentError forcedensity_bc!((p,t) -> p[1] * t, body, :b, :x)

    data = rand(2, 2)
    Peridynamics.forcedensity_databc!(body, data, :a, [:y, :z])
    @test_throws ArgumentError forcedensity_bc!(t -> t, body, :a, :y)
    @test_throws ArgumentError forcedensity_bc!(t -> t, body, :b, :y)
    @test_throws ArgumentError forcedensity_bc!(t -> t, body, :a, :z)
    @test_throws ArgumentError forcedensity_bc!(t -> t, body, :b, :z)

    # unknown symbols
    data = rand(1, 2)
    @test_throws ArgumentError Peridynamics.velocity_databc!(body, data, :a, [:k])
    @test_throws ArgumentError Peridynamics.velocity_databc!(body, data, :a, [4])
    @test_throws MethodError Peridynamics.velocity_databc!(body, data, :a, 3)
    data = rand(2, 2)
    @test_throws MethodError Peridynamics.velocity_databc!(body, data, :a, [:x, 2])

    # too many dimensions
    data = rand(4, 2)
    @test_throws ArgumentError Peridynamics.velocity_databc!(body, data, :a, [1, 2, 3, 4])
end
