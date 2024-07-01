@testitem "show Job" begin
    io = IOBuffer()

    b1 = Body(BBMaterial(), rand(3,10), rand(10))
    material!(b1, horizon=1, E=1, rho=1, Gc=1)
    velocity_ic!(b1, :all_points, :x, 1.0)
    b2 = Body(OSBMaterial(), rand(3,5), rand(5))
    material!(b2, horizon=1, E=1, nu=0.25, rho=1, Gc=1)
    ms = MultibodySetup(:a => b1, :b => b2)
    job = Job(ms, VelocityVerlet(steps=1))

    show(IOContext(io, :compact=>true), MIME("text/plain"), job)
    msg = String(take!(io))
    @test contains(msg, "15-point multibody Job with VelocityVerlet solver")

    show(IOContext(io, :compact=>false), MIME("text/plain"), job)
    msg = String(take!(io))
    @test contains(msg, "15-point MultibodySetup")
    @test contains(msg, "VelocityVerlet(n_steps=1, safety_factor=0.7)")
end

@testitem "Job pre submission checks" begin
    b1 = Body(BBMaterial(), rand(3,10), rand(10))
    vv = VelocityVerlet(steps=1)
    @test_throws ErrorException Job(b1, vv)

    point_set!(b1, :a, 1:2)
    material!(b1, :a; horizon=1, E=1, rho=1, Gc=1)
    @test_throws ErrorException Job(b1, vv)

    material!(b1, :all_points; horizon=1, E=1, rho=1, Gc=1)
    velocity_bc!(t -> 0, b1, :all_points, 1)
    job = Job(b1, vv)
    @test job.spatial_setup isa Body{<:BBMaterial}

    b2 = Body(BBMaterial(), rand(3,10), rand(10))
    b3 = Body(OSBMaterial(), rand(3,5), rand(5))
    ms = MultibodySetup(:b2 => b2, :b3 => b3)
    @test_throws ErrorException Job(ms, vv)

    material!(b2, horizon=1, E=1, rho=1, Gc=1)
    material!(b3, horizon=1, E=1, nu=0.25, rho=1, Gc=1)
    @test_throws ErrorException Job(ms, vv)

    velocity_ic!(b2, :all_points, 1, 1.0)
    job = Job(ms, vv)
    @test job.spatial_setup isa MultibodySetup
end
