@testitem "boundary conditions" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    point_set!(body, :a, 1:2)
    velocity_bc!(t -> t, body, :a, :x)
    forcedensity_bc!(t -> t, body, :a, :x)
    ts = VelocityVerlet(steps=10)
    pd = Peridynamics.PointDecomposition(body, 2)
    ps = Peridynamics.get_param_spec(body)
    chunk = Peridynamics.BodyChunk(body, ts, pd, 1, ps)
    Peridynamics.apply_bcs!(chunk, 1.0)

    (; velocity_half, b_ext) = chunk.storage
    @test chunk.storage.velocity_half ≈ [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test chunk.storage.b_ext ≈ [1.0 1.0; 0.0 0.0; 0.0 0.0]
end

@testitem "test all conditions" begin
    position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=1.5, rho=1, E=210e9, Gc=1.0)
    point_set!(body, :a, 1:2)
    velocity_bc!(t -> t, body, :a, :x)
    forcedensity_bc!(t -> t, body, :a, :x)
    velocity_ic!(body, :a, :x, 1.234)
    velocity_ic!(p -> 3.456, body, :a, :y)
    @test_throws ArgumentError velocity_ic!(t -> 1.0, body, :a, :z)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; velocity_half, velocity, b_ext) = chunk.storage

    @test velocity ≈ [1.234 1.234; 3.456 3.456; 0.0 0.0]

    Peridynamics.apply_bcs!(chunk, 2.345)

    @test velocity_half ≈ [2.345 2.345; 0.0 0.0; 0.0 0.0]
    @test b_ext ≈ [2.345 2.345; 0.0 0.0; 0.0 0.0]
end

@testitem "initial conditions" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    point_set!(body, :a, 1:2)
    velocity_bc!(t -> t, body, :a, :x)
    forcedensity_bc!(t -> t, body, :a, :x)
    ts = VelocityVerlet(steps=10)
    pd = Peridynamics.PointDecomposition(body, 2)
    ps = Peridynamics.get_param_spec(body)
    chunk = Peridynamics.BodyChunk(body, ts, pd, 1, ps)
    Peridynamics.apply_bcs!(chunk, 1.0)

    (; velocity_half, b_ext) = chunk.storage
    @test chunk.storage.velocity_half ≈ [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test chunk.storage.b_ext ≈ [1.0 1.0; 0.0 0.0; 0.0 0.0]
end
