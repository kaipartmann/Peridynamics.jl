@testitem "MultibodySetup" begin
    b1 = Body(BBMaterial(), rand(3,10), rand(10))
    b2 = Body(BBMaterial(), rand(3,10), rand(10))
    b3 = Body(OSBMaterial(), rand(3,10), rand(10))

    ms = MultibodySetup(:a => b1, :b => b2)
    @test ms.bodies == (b1, b2) ||  ms.bodies == (b2, b1)
    @test ms.body_names == [:a, :b] || ms.body_names == [:b, :a]
    @test ms.body_idxs == Dict(:a => 1, :b => 2) || ms.body_idxs == Dict(:a => 2, :b => 1)
    @test isempty(ms.srf_contacts)

    ms = MultibodySetup(Dict(:a => b1, :b => b2))
    @test ms.bodies == (b1, b2) ||  ms.bodies == (b2, b1)
    @test ms.body_names == [:a, :b] || ms.body_names == [:b, :a]
    @test ms.body_idxs == Dict(:a => 1, :b => 2) || ms.body_idxs == Dict(:a => 2, :b => 1)
    @test isempty(ms.srf_contacts)

    ms = MultibodySetup(:a => b1, :b => b2, :c => b3)
    @test ms.bodies[ms.body_idxs[:a]] == b1
    @test ms.bodies[ms.body_idxs[:b]] == b2
    @test ms.bodies[ms.body_idxs[:c]] == b3
    for name in (:a, :b, :c)
        idx = ms.body_idxs[name]
        @test Peridynamics.get_body(ms, name) == ms.bodies[idx]
        @test Peridynamics.get_body_name(ms, idx) == string(name)
    end

    @test_throws ArgumentError MultibodySetup(:a => b1)
    @test_throws ArgumentError MultibodySetup(Dict(:a => b1))
end

@testitem "contact!" begin
    using Peridynamics.PointNeighbors

    b1 = Body(BBMaterial(), rand(3,10), rand(10))
    b2 = Body(BBMaterial(), rand(3,10), rand(10))
    ms = MultibodySetup(:body1 => b1, :body2 => b2)

    nhs = GridNeighborhoodSearch{3}(search_radius=1, n_points=10)

    contact!(ms, :body1, :body2; radius=1)
    srfc = ms.srf_contacts[1]
    @test srfc.body_id_a === :body1
    @test srfc.body_id_b === :body2
    @test srfc.radius == 1.0
    @test srfc.penalty_factor == 1e12
    @test srfc.nhs isa GridNeighborhoodSearch{3}
    empty!(ms.srf_contacts)

    contact!(ms, :body1, :body2; radius=2, penalty_factor=1.5e12)
    srfc = ms.srf_contacts[1]
    @test srfc.body_id_a === :body1
    @test srfc.body_id_b === :body2
    @test srfc.radius == 2.0
    @test srfc.penalty_factor == 1.5e12
    @test srfc.nhs isa GridNeighborhoodSearch{3}
    empty!(ms.srf_contacts)

    @test_throws UndefKeywordError(:radius) contact!(ms, :body1, :body2)

    @test_throws ArgumentError contact!(ms, :body1, :body2; radius=0)

    @test_throws ArgumentError contact!(ms, :body1, :body2; radius=1, penalty_factor=0)

    @test_throws ArgumentError contact!(ms, :body1, :b; radius=1)
end

@testitem "show MultibodySetup" begin
    io = IOBuffer()

    b1 = Body(BBMaterial(), rand(3,10), rand(10))
    material!(b1, horizon=1, E=1, rho=1, Gc=1)
    velocity_ic!(b1, :all_points, :x, 1.0)
    b2 = Body(OSBMaterial(), rand(3,5), rand(5))
    material!(b2, horizon=1, E=1, nu=0.25, rho=1, Gc=1)
    ms = MultibodySetup(:a => b1, :b => b2)

    show(IOContext(io, :compact=>true), MIME("text/plain"), ms)
    msg = String(take!(io))
    @test contains(msg, "15-point MultibodySetup")

    show(IOContext(io, :compact=>false), MIME("text/plain"), ms)
    msg = String(take!(io))
    @test contains(msg, "15-point MultibodySetup")
    @test contains(msg, "10-point Body")
    @test contains(msg, "BBMaterial")
    @test contains(msg, "with name `a`")
    @test contains(msg, "5-point Body")
    @test contains(msg, "OSBMaterial")
    @test contains(msg, "with name `b`")
end
