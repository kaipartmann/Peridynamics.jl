@testitem "MultibodySetup" begin
    b1 = Body(BBMaterial(), rand(3,10), rand(10))
    b2 = Body(BBMaterial(), rand(3,10), rand(10))
    b3 = Body(CKIMaterial(), rand(3,10), rand(10))

    ms = MultibodySetup(:body1 => b1, :body2 => b2)
    @test ms isa MultibodySetup{Peridynamics.BBMaterial{Peridynamics.NoCorrection},Peridynamics.BBPointParameters}
    @test Peridynamics.material_type(ms) == Peridynamics.BBMaterial{Peridynamics.NoCorrection}
    @test length(ms.bodies) == 2
    @test isempty(ms.srf_contacts)

    ms = MultibodySetup(Dict(:body1 => b1, :body2 => b2))
    @test ms isa MultibodySetup{Peridynamics.BBMaterial{Peridynamics.NoCorrection},Peridynamics.BBPointParameters}
    @test Peridynamics.material_type(ms) == Peridynamics.BBMaterial{Peridynamics.NoCorrection}
    @test length(ms.bodies) == 2
    @test isempty(ms.srf_contacts)

    @test_throws ArgumentError MultibodySetup(:body1 => b1, :body3 => b3)
end

@testitem "contact!" begin
    using Peridynamics.PointNeighbors

    b1 = Body(BBMaterial(), rand(3,10), rand(10))
    b2 = Body(BBMaterial(), rand(3,10), rand(10))
    ms = MultibodySetup(:body1 => b1, :body2 => b2)

    nhs = GridNeighborhoodSearch{3}(1, 10)

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
end
