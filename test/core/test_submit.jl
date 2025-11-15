@testitem "get_submit_options" begin
    o = Dict{Symbol,Any}(:quiet => true)
    quiet = Peridynamics.get_submit_options(o)
    @test quiet == true

    o = Dict{Symbol,Any}()
    quiet = Peridynamics.get_submit_options(o)
    @test quiet == false
end

@testitem "submit threads error handling" begin
    ref_position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    body = Body(BBMaterial(), ref_position, volume)
    material!(body, horizon=1.5, rho=8000, E=210e9)
    Δt = 1e-7
    velocity_bc!(t -> t < 1.9Δt ? 1.0 : 1e40, body, :all_points, :x)
    ts = VelocityVerlet(steps=5, stepsize=Δt)
    job = Job(body, ts; path= mktempdir())
    @test_throws CompositeException Peridynamics.submit_threads(job, 2)
    logfile_contents = read(joinpath(job.options.logfile), String)
    @test occursin("NaN values detected", logfile_contents)
    @test occursin("Stacktrace", logfile_contents)
end

@testitem "submit MPI error handling" begin
    #-- simulation should fail at step 2 due to NaN values --#
    pos = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    vol = [1.0, 1.0]
    body = Body(BBMaterial(), pos, vol)
    material!(body, horizon=1.5, rho=8000, E=210e9)
    Δt = 1e-7
    velocity_bc!(t -> t < 1.9Δt ? 1.0 : 1e40, body, :all_points, :x)
    ts = VelocityVerlet(steps=5, stepsize=Δt)
    job = Job(body, ts; path= mktempdir())
    err = Peridynamics.NaNError(2Δt, 2) # extected error
    @test_throws err Peridynamics.submit_mpi(job) # error is thrown!
    logfile_contents = read(joinpath(job.options.logfile), String)
    @test occursin("NaN values detected", logfile_contents)
    @test occursin("Stacktrace", logfile_contents)

    #-- simulation should fail at step 2 due to some weird error --#
    pos = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    vol = [1.0, 1.0]
    body = Body(BBMaterial(), pos, vol)
    material!(body, horizon=1.5, rho=8000, E=210e9)
    Δt = 1e-7
    err_msg = "some weird error occurred!\n"
    err = ErrorException(err_msg)
    velocity_bc!(body, :all_points, :x) do t
        t < 1.9Δt ? 1.0 : throw(err)
    end
    ts = VelocityVerlet(steps=5, stepsize=Δt)
    job = Job(body, ts; path= mktempdir())
    @test_throws err Peridynamics.submit_mpi(job) # the same error should be thrown here!
    logfile_contents = read(joinpath(job.options.logfile), String)
    @test occursin(err_msg, logfile_contents)
end

@testitem "submit MPI NaNError with multiple ranks" tags=[:mpi] begin
    path = mktempdir()
    Δt = 1e-7
    mpi_cmd = """
    using Peridynamics, Test
    #-- simulation should fail at step 2 due to NaN values --#
    pos = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    vol = [1.0, 1.0]
    body = Body(BBMaterial(), pos, vol)
    material!(body, horizon=1.5, rho=8000, E=210e9)
    Δt = $(Δt)
    velocity_bc!(t -> t < 1.9Δt ? 1.0 : 1e40, body, :all_points, :x)
    ts = VelocityVerlet(steps=5, stepsize=Δt)
    job = Job(body, ts; path="$(path)")
    err = Peridynamics.NaNError(2Δt, 2) # extected error
    @test_throws err submit(job) # extected error is thrown!
    """
    mpiexec = Peridynamics.MPI.mpiexec()
    jlcmd = Base.julia_cmd()
    pdir = pkgdir(Peridynamics)
    cmd = `$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd)`
    @test success(cmd) # does not print anything
    # for debugging use the run command:
    # run(cmd)

    logfile_contents = read(joinpath(path, "logfile.log"), String)
    msg_correct = "NaN values detected in simulation data!\n"
    msg_correct *= "  time:    $(2Δt)\n  step:    2\n"
    @test occursin(msg_correct, logfile_contents)
end

@testitem "submit MPI Abort with multiple ranks" tags=[:mpi] begin
    path = mktempdir()
    Δt = 1e-7
    mpi_cmd = """
    using Peridynamics, Test
    #-- simulation should fail at step 2 due to NaN values --#
    pos = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    vol = [1.0, 1.0]
    body = Body(BBMaterial(), pos, vol)
    material!(body, horizon=1.5, rho=8000, E=210e9)
    Δt = $(Δt)
    velocity_bc!(body, :all_points, :x) do t
        t < 1.9Δt ? 1.0 : error("some weird error occurred!\n")
    end
    ts = VelocityVerlet(steps=5, stepsize=Δt)
    job = Job(body, ts; path="$(path)")
    @test_throws ErrorException submit(job) # The whole thing is aborted!
    """
    mpiexec = Peridynamics.MPI.mpiexec()
    jlcmd = Base.julia_cmd()
    pdir = pkgdir(Peridynamics)
    cmd = `$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd)`
    @test !success(cmd) # does not print anything
    # for debugging use the run command:
    # run(cmd)

    logfile_contents = read(joinpath(path, "logfile.log"), String)
    @test occursin("some weird error occurred!", logfile_contents)
end
