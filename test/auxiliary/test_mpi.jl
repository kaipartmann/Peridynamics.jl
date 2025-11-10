@testitem "force_mpi_run!" begin
    Peridynamics.force_mpi_run!()
    @test Peridynamics.mpi_run() == true
    Peridynamics.init_mpi() # will do nothing because MPI_RUN was forced
    @test Peridynamics.mpi_run() == true
end

@testitem "force_threads_run!" begin
    Peridynamics.force_threads_run!()
    @test Peridynamics.mpi_run() == false
    Peridynamics.init_mpi() # will do nothing because MPI_RUN was forced
    @test Peridynamics.mpi_run() == false
end

@testitem "@mpiroot error case" begin
    @test_throws LoadError eval(:(@mpiroot :wrongkey println("Hi!")))
end

@testitem "mpi_barrier" begin
    # Test that mpi_barrier returns nothing
    @test mpi_barrier() === nothing
end

@testitem "MPI timers" tags=[:mpi] begin
    path = mktempdir()
    mpi_cmd = """
    using Peridynamics
    disable_mpi_timers!() # disable first
    enable_mpi_timers!()
    path = "$(path)"
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = CMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, nu=0.25, Gc=1)
    velocity_ic!(body, :all_points, :x, 1.0)
    ts = VelocityVerlet(steps=10)
    job = Job(body, ts; path, freq=10)
    submit(job)
    """
    mpiexec = Peridynamics.MPI.mpiexec()
    jlcmd = Base.julia_cmd()
    pdir = pkgdir(Peridynamics)
    cmd = `$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd)`
    @test success(cmd) # does not print anything
    # for debugging use the run command:
    # run(cmd)

    # test if the timer logfiles exists
    created_files = filter(x -> endswith(x, ".log"), readdir(path))
    @test any(x -> x == "timers_rank_0.log", created_files)
    @test any(x -> x == "timers_rank_1.log", created_files)
end
