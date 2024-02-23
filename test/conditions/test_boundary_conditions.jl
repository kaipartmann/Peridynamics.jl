@testitem "apply_bc!" begin
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
    b = Peridynamics.BodyChunk(body, ts, pd, 1)
    Peridynamics.apply_bcs!(b, 1.0)

    @test b.store.velocity_half ≈ [1.0 1.0; 0.0 0.0; 0.0 0.0]
    @test b.store.b_ext ≈ [1.0 1.0; 0.0 0.0; 0.0 0.0]

end
