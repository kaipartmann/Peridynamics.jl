@testitem "MultibodySetup" begin
    b1 = Body(BBMaterial(), rand(3,10), rand(10))
    b2 = Body(BBMaterial(), rand(3,10), rand(10))
    b3 = Body(CKIMaterial(), rand(3,10), rand(10))

    ms = MultibodySetup(:body1 => b1, :body2 => b2)
    @test ms isa MultibodySetup{Peridynamics.BBMaterial,Peridynamics.BBPointParameters}
    @test Peridynamics.material_type(ms) == Peridynamics.BBMaterial
    @test length(ms.bodies) == 2
    @test isempty(ms.contacts)

    ms = MultibodySetup(Dict(:body1 => b1, :body2 => b2))
    @test ms isa MultibodySetup{Peridynamics.BBMaterial,Peridynamics.BBPointParameters}
    @test Peridynamics.material_type(ms) == Peridynamics.BBMaterial
    @test length(ms.bodies) == 2
    @test isempty(ms.contacts)

    @test_throws ArgumentError MultibodySetup(:body1 => b1, :body3 => b3)
end

@testitem "contact!" begin
    b1 = Body(BBMaterial(), rand(3,10), rand(10))
    b2 = Body(BBMaterial(), rand(3,10), rand(10))
    ms = MultibodySetup(:body1 => b1, :body2 => b2)

    contact!(ms, :body1, :body2; radius=1)
    @test ms.contacts[1] == Peridynamics.Contact(:body1, :body2, 1.0, 1e12)
    empty!(ms.contacts)

    contact!(ms, :body1, :body2; radius=2, sc=1.5e12)
    @test ms.contacts[1] == Peridynamics.Contact(:body1, :body2, 2.0, 1.5e12)
    empty!(ms.contacts)

    @test_throws UndefKeywordError(:radius) contact!(ms, :body1, :body2)

    @test_throws ArgumentError contact!(ms, :body1, :body2; radius=0)

    @test_throws ArgumentError contact!(ms, :body1, :body2; radius=1, sc=0)
end
