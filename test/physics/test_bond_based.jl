@testitem "get_storage" begin
    position, volume = uniform_box(1,1,1,0.25)
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    ts = VelocityVerlet(steps=10)
    pd = Peridynamics.PointDecomposition(body, 1)
    db, ch = Peridynamics.get_system(body, pd, 1)
    s = Peridynamics.get_storage(body.mat, ts, db, ch)

    @test s isa Peridynamics.BBVerletStorage
    @test s.position == position
    @test s.displacement == zeros(3, 64)
    @test s.velocity == zeros(3, 64)
    @test s.velocity_half == zeros(3, 64)
    @test s.acceleration == zeros(3, 64)
    @test s.b_int == zeros(3, 64)
    @test s.b_ext == zeros(3, 64)
    @test s.damage == zeros(64)
    @test s.bond_active == ones(Bool, 4032)
    @test s.n_active_bonds == fill(63, 64)
end
