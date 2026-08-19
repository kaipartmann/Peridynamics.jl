@testitem "show Job" begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

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

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "Job pre submission checks" begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

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

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "semidiscretize Job" begin
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    body = Body(BBMaterial(), position, ones(2))
    material!(body; horizon=1.5, E=1.0, rho=1.0, Gc=1.0)
    velocity_ic!(body, :all_points, :x, 1.0)
    job = Job(body, VelocityVerlet(steps=2, stepsize=0.1))

    ode = semidiscretize(job; n_chunks=1)
    @test ode isa Peridynamics.SciMLBase.ODEProblem
    @test ode.tspan == (0.0, 0.2)
    @test ode.p.job === job
    @test ode.p.data_handler isa Peridynamics.AbstractThreadsBodyDataHandler
    @test length(ode.u0.x[1]) == 6
    @test length(ode.u0.x[2]) == 6

    du = similar(ode.u0)
    ode.f(du, ode.u0, ode.p, 0.0)
    @test du.x[1] ≈ zeros(6)
    @test du.x[2] ≈ [1.0, 0.0, 0.0, 1.0, 0.0, 0.0]

    Peridynamics.MPI_RUN[] = mpi_run_current_value
end
