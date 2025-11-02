@testitem "Study constructor with valid setups" begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    function create_job(setup::NamedTuple)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=1, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, setup.velocity)
        vv = VelocityVerlet(steps=setup.n_steps)
        job = Job(body, vv)
        return job
    end

    setups = [
        (; n_steps=10, velocity=1.0),
        (; n_steps=20, velocity=2.0),
        (; n_steps=30, velocity=3.0),
    ]

    study = Study(create_job, setups)

    @test study isa Study
    @test length(study.jobs) == 3
    @test length(study.setups) == 3
    @test length(study.submission_status) == 3
    @test length(study.postproc_status) == 3
    @test length(study.results) == 3
    @test all(.!study.submission_status)
    @test all(.!study.postproc_status)
    @test all(isempty.(study.results))

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "Study constructor with empty setups" begin
    function create_job(setup::NamedTuple)
        return Job(Body(BBMaterial(), rand(3, 10), rand(10)), VelocityVerlet(steps=1))
    end

    setups = NamedTuple[]

    @test_throws ArgumentError Study(create_job, setups)
end

@testitem "Study constructor with inconsistent setup fields" begin
    function create_job(setup::NamedTuple)
        return Job(Body(BBMaterial(), rand(3, 10), rand(10)), VelocityVerlet(steps=1))
    end

    setups = [
        (; n_steps=10, velocity=1.0),
        (; n_steps=20, force=2.0),  # Different field name
    ]

    @test_throws ArgumentError Study(create_job, setups)
end

@testitem "Study constructor with create_job that doesn't return Job" begin
    function create_job(setup::NamedTuple)
        return "not a job"
    end

    setups = [(; n_steps=10,)]

    @test_throws ArgumentError Study(create_job, setups)
end

@testitem "Study constructor with create_job that errors" begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    function create_job(setup::NamedTuple)
        if setup.n_steps < 0
            error("Invalid n_steps")
        end
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=1, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        return Job(body, VelocityVerlet(steps=setup.n_steps))
    end

    setups = [
        (; n_steps=10,),
        (; n_steps=-5,),
    ]

    @test_throws ArgumentError Study(create_job, setups)

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "Study show methods" begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    function create_job(setup::NamedTuple)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=1, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        return Job(body, VelocityVerlet(steps=setup.n_steps))
    end

    setups = [
        (; n_steps=10,),
        (; n_steps=20,),
    ]

    study = Study(create_job, setups)

    io = IOBuffer()
    show(IOContext(io, :compact=>true), study)
    msg = String(take!(io))
    @test contains(msg, "Study with 2 simulations")
    @test contains(msg, "0 submitted")
    @test contains(msg, "0 post-processed")

    show(IOContext(io, :compact=>false), MIME("text/plain"), study)
    msg = String(take!(io))
    @test contains(msg, "Study:")
    @test contains(msg, "Number of simulations:")
    @test contains(msg, "Setup parameters:")

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "submit! with all successful jobs" tags=[:skipci] begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, setup.velocity)
            vv = VelocityVerlet(steps=setup.n_steps)
            path = joinpath(tmpdir, "sim_$(setup.n_steps)")
            job = Job(body, vv; path=path, freq=5)
            return job
        end

        setups = [
            (; n_steps=5, velocity=1.0),
            (; n_steps=10, velocity=2.0),
        ]

        study = Study(create_job, setups)

        @test all(.!study.submission_status)

        submit!(study; quiet=true)

        @test all(study.submission_status)
        @test count(study.submission_status) == 2
    end

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "submit! with some failing jobs" tags=[:skipci] begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, setup.velocity)
            vv = VelocityVerlet(steps=setup.n_steps)

            # Create a job with invalid path for the second simulation
            # This will fail during submission when trying to create directories
            if setup.fail
                path = "/invalid/path/that/cannot/be/created/sim_$(setup.n_steps)"
            else
                path = joinpath(tmpdir, "sim_$(setup.n_steps)")
            end

            job = Job(body, vv; path=path, freq=5)
            return job
        end

        setups = [
            (; n_steps=5, velocity=1.0, fail=false),
            (; n_steps=10, velocity=2.0, fail=true),
            (; n_steps=15, velocity=3.0, fail=false),
        ]

        study = Study(create_job, setups)

        # This should not throw, but should mark the second job as failed
        submit!(study; quiet=true)

        @test study.submission_status[1] == true
        @test study.submission_status[2] == false
        @test study.submission_status[3] == true
        @test count(study.submission_status) == 2
    end

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "postproc! without prior submission" tags=[:skipci] begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            vv = VelocityVerlet(steps=5)
            path = joinpath(tmpdir, "sim_$(setup.id)")
            job = Job(body, vv; path=path, freq=5)
            return job
        end

        setups = [(; id=1,)]
        study = Study(create_job, setups)

        function proc_func(r0, r, id)
            return nothing
        end

        # Should warn and return nothing
        postproc!(proc_func, study)
        @test isempty(study.results[1])
    end

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "postproc! with successful submission returning nothing" tags=[:skipci] begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            vv = VelocityVerlet(steps=5)
            path = joinpath(tmpdir, "sim_$(setup.id)")
            job = Job(body, vv; path=path, freq=5)
            return job
        end

        setups = [(; id=1,)]
        study = Study(create_job, setups)
        submit!(study; quiet=true)

        function proc_func(r0, r, id)
            # Do some processing but return nothing
            return nothing
        end

        postproc!(proc_func, study; serial=true)
        @test study.postproc_status[1] == true
        @test isempty(study.results[1])
    end

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "postproc! with successful submission returning NamedTuples" tags=[:skipci] begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            vv = VelocityVerlet(steps=10)
            path = joinpath(tmpdir, "sim_$(setup.id)")
            job = Job(body, vv; path=path, freq=5)
            return job
        end

        setups = [
            (; id=1,),
            (; id=2,),
        ]
        study = Study(create_job, setups)
        submit!(study; quiet=true)

        function proc_func(r0, r, id)
            max_disp = maximum(abs.(r[:displacement]))
            return (; time_step=id, max_displacement=max_disp)
        end

        postproc!(proc_func, study; serial=true)

        @test all(study.postproc_status)
        @test !isempty(study.results[1])
        @test !isempty(study.results[2])
        @test all(r -> r isa NamedTuple, study.results[1])
        @test all(r -> r isa NamedTuple, study.results[2])
        @test all(r -> haskey(r, :time_step) && haskey(r, :max_displacement), study.results[1])
    end

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "postproc! with some jobs having errors" tags=[:skipci] begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            vv = VelocityVerlet(steps=10)
            path = joinpath(tmpdir, "sim_$(setup.id)")
            job = Job(body, vv; path=path, freq=5)
            return job
        end

        setups = [
            (; id=1,),
            (; id=2,),
        ]
        study = Study(create_job, setups)
        submit!(study; quiet=true)

        # Delete the vtk directory for the first simulation to cause an error
        rm(joinpath(tmpdir, "sim_1", "vtk"); recursive=true)

        function proc_func(r0, r, id)
            max_disp = maximum(abs.(r[:displacement]))
            return (; time_step=id, max_displacement=max_disp)
        end

        postproc!(proc_func, study; serial=true)

        @test study.postproc_status[1] == false
        @test study.postproc_status[2] == true
        @test count(study.postproc_status) == 1
        @test isempty(study.results[1])
        @test !isempty(study.results[2])
    end

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "postproc! with invalid processing function" tags=[:skipci] begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            vv = VelocityVerlet(steps=5)
            path = joinpath(tmpdir, "sim_$(setup.id)")
            job = Job(body, vv; path=path, freq=5)
            return job
        end

        setups = [(; id=1,)]
        study = Study(create_job, setups)
        submit!(study; quiet=true)

        # Function with wrong number of arguments
        function bad_proc_func(r0, r)
            return nothing
        end

        @test_throws ArgumentError postproc!(bad_proc_func, study)
    end

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "Study with MultibodySetup" tags=[:skipci] begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple)
            b1 = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(b1, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(b1, :all_points, :x, setup.velocity)

            b2 = Body(OSBMaterial(), rand(3, 5), rand(5))
            material!(b2, horizon=1, E=1, nu=0.25, rho=1, Gc=1)

            ms = MultibodySetup(:body1 => b1, :body2 => b2)
            vv = VelocityVerlet(steps=setup.n_steps)
            job = Job(ms, vv)
            return job
        end

        setups = [
            (; id=1, n_steps=5, velocity=1.0),
        ]

        study = Study(create_job, setups)
        @test length(study.jobs) == 1
        @test study.jobs[1].spatial_setup isa MultibodySetup
    end

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end
