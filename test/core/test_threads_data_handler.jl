@testitem "ThreadsBodyDataHandler BBMaterial VelocityVerlet" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    point_set!(body, :a, 1:2)
    point_set!(body, :b, 3:4)
    velocity_ic!(body, :a, :x, 1.0)
    velocity_bc!(t->t, body, :a, :x)
    forcedensity_bc!(t->t, body, :a, :x)
    precrack!(body, :a, :b)
    ts = VelocityVerlet(steps=10)
    tdh = Peridynamics.ThreadsBodyDataHandler(body, ts, 2)

    @test tdh.n_chunks == 2
    @test length(tdh.chunks) == 2
    @test length(tdh.lth_exs[1]) == 1
    @test length(tdh.lth_exs[2]) == 1
    @test length(tdh.htl_exs[1]) == 1
    @test length(tdh.htl_exs[2]) == 1

    #TODO
end
