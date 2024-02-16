@testitem "init_storage" begin
    body = Body(BBMaterial(), rand(3,10), rand(10))
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    ts = VelocityVerlet(steps=10)
    db, ch = Peridynamics.init_discretization(body, 1:10)
    s = Peridynamics.init_storage(body, ts, db, ch)
end
