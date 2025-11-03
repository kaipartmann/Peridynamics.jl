@testitem "Study constructor with valid setups" begin
    # for now, this test only works with multithreading!
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple, root::String)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=setup.E, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            vv = VelocityVerlet(steps=setup.n_steps)
            path = joinpath(root, "sim_$(setup.E)_$(setup.n_steps)")
            job = Job(body, vv; path=path, freq=5)
            return job
        end

        setups = [
            (; E=1.0, n_steps=5),
            (; E=2.0, n_steps=10),
            (; E=3.0, n_steps=15),
        ]

        study = Study(create_job, setups; root=joinpath(tmpdir, "study"))

        @test study isa Study
        @test length(study.jobs) == 3
        @test length(study.setups) == 3
        @test length(study.sim_success) == 3
        @test length(study.jobpaths) == 3
        @test all(.!study.sim_success)
        @test study.root == joinpath(tmpdir, "study")
        @test study.logfile == joinpath(tmpdir, "study", "study_log.log")
        @test all(contains.(study.jobpaths, "sim_"))
    end

    # reset to the value as before
    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "Study constructor with empty setups" begin
    mktempdir() do tmpdir
        function create_job(setup::NamedTuple, root::String)
            return Job(Body(BBMaterial(), rand(3, 10), rand(10)), VelocityVerlet(steps=1))
        end

        setups = NamedTuple[]

        @test_throws ArgumentError Study(create_job, setups; root=tmpdir)
    end
end

@testitem "Study constructor with inconsistent setup fields" begin
    mktempdir() do tmpdir
        function create_job(setup::NamedTuple, root::String)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            return Job(body, VelocityVerlet(steps=1))
        end

        setups = [
            (; n_steps=10, velocity=1.0),
            (; n_steps=20, force=2.0),  # Different field name
        ]

        @test_throws ArgumentError Study(create_job, setups; root=tmpdir)
    end
end

@testitem "Study constructor with jobcreator that errors" begin
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple, root::String)
            if setup.n_steps < 0
                error("Invalid n_steps")
            end
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            path = joinpath(root, "sim_$(setup.n_steps)")
            return Job(body, VelocityVerlet(steps=setup.n_steps); path=path)
        end

        setups = [
            (; n_steps=10,),
            (; n_steps=-5,),
        ]

        @test_throws ErrorException Study(create_job, setups; root=tmpdir)
    end

    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "Study constructor with duplicate job paths" begin
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple, root::String)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            # All jobs use the same path - should error!
            path = joinpath(root, "sim")
            return Job(body, VelocityVerlet(steps=setup.n_steps); path=path)
        end

        setups = [
            (; n_steps=10,),
            (; n_steps=20,),
        ]

        @test_throws ArgumentError Study(create_job, setups; root=tmpdir)
    end

    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "Study show method" begin
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple, root::String)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            path = joinpath(root, "sim_$(setup.n_steps)")
            return Job(body, VelocityVerlet(steps=setup.n_steps); path=path)
        end

        setups = [
            (; n_steps=10,),
            (; n_steps=20,),
        ]

        study = Study(create_job, setups; root=tmpdir)

        io = IOBuffer()
        show(io, MIME("text/plain"), study)
        msg = String(take!(io))

        @test contains(msg, "Study with 2 jobs")
        @test contains(msg, "✗")  # Not yet submitted
        @test contains(msg, "sim_10")
        @test contains(msg, "sim_20")
    end

    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "submit! with all successful jobs" tags=[:skipci] begin
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple, root::String)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, setup.velocity)
            vv = VelocityVerlet(steps=setup.n_steps)
            path = joinpath(root, "sim_$(setup.n_steps)")
            job = Job(body, vv; path=path, freq=5)
            return job
        end

        setups = [
            (; n_steps=5, velocity=1.0),
            (; n_steps=10, velocity=2.0),
        ]

        study = Study(create_job, setups; root=joinpath(tmpdir, "study"))

        @test all(.!study.sim_success)

        submit!(study; quiet=true)

        @test all(study.sim_success)
        @test count(study.sim_success) == 2

        # Check that logfile was created
        @test isfile(study.logfile)
        logcontent = read(study.logfile, String)
        @test contains(logcontent, "SIMULATION STUDY LOGFILE")
        @test contains(logcontent, "sim_5")
        @test contains(logcontent, "sim_10")
        @test contains(logcontent, "completed ✓")
        @test contains(logcontent, "n_steps: 5")
        @test contains(logcontent, "velocity: 1.0")
    end

    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "submit! with some failing jobs" tags=[:skipci] begin
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple, root::String)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            vv = VelocityVerlet(steps=setup.n_steps)

            # Create invalid path for second job
            if setup.fail
                path = "/invalid/path/that/cannot/be/created/sim_$(setup.n_steps)"
            else
                path = joinpath(root, "sim_$(setup.n_steps)")
            end

            job = Job(body, vv; path=path, freq=5)
            return job
        end

        setups = [
            (; n_steps=5, fail=false),
            (; n_steps=10, fail=true),
            (; n_steps=15, fail=false),
        ]

        study = Study(create_job, setups; root=joinpath(tmpdir, "study"))

        submit!(study; quiet=true)

        @test study.sim_success[1] == true
        @test study.sim_success[2] == false
        @test study.sim_success[3] == true
        @test count(study.sim_success) == 2

        # Check logfile contains failure info
        logcontent = read(study.logfile, String)
        @test contains(logcontent, "completed ✓")
        @test contains(logcontent, "failed ✗")
    end

    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "submit! creates study directory structure" tags=[:skipci] begin
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple, root::String)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            vv = VelocityVerlet(steps=setup.n_steps)
            path = joinpath(root, "sim_$(setup.n_steps)")
            job = Job(body, vv; path=path, freq=5)
            return job
        end

        setups = [(; n_steps=5,)]
        studyroot = joinpath(tmpdir, "my_study")

        study = Study(create_job, setups; root=studyroot)

        @test !isdir(studyroot)

        submit!(study; quiet=true)

        @test isdir(studyroot)
        @test isfile(study.logfile)
        @test isdir(study.jobpaths[1])
        @test isfile(joinpath(study.jobpaths[1], "logfile.log"))
    end

    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "submit! with quiet keyword" tags=[:skipci] begin
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple, root::String)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            vv = VelocityVerlet(steps=setup.n_steps)
            path = joinpath(root, "sim_$(setup.n_steps)")
            job = Job(body, vv; path=path, freq=5)
            return job
        end

        setups = [(; n_steps=5,)]
        study = Study(create_job, setups; root=joinpath(tmpdir, "study"))

        # Should not error with quiet=true
        submit!(study; quiet=true)
        @test study.sim_success[1] == true

        # Should not error with quiet=false (default)
        study2 = Study(create_job, setups; root=joinpath(tmpdir, "study2"))
        submit!(study2)
        @test study2.sim_success[1] == true
    end

    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "Study with MultibodySetup" tags=[:skipci] begin
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple, root::String)
            b1 = Body(BBMaterial(), rand(3, 50), rand(50))
            material!(b1, horizon=1, E=1, rho=1, Gc=1)
            velocity_ic!(b1, :all_points, :x, setup.velocity)

            b2 = Body(OSBMaterial(), rand(3, 50), rand(50))
            material!(b2, horizon=1, E=1, nu=0.25, rho=1, Gc=1)

            ms = MultibodySetup(:body1 => b1, :body2 => b2)
            vv = VelocityVerlet(steps=setup.n_steps)
            path = joinpath(root, "sim_$(setup.id)")
            job = Job(ms, vv; path=path, freq=5)
            return job
        end

        setups = [(; id=1, n_steps=5, velocity=1.0)]
        study = Study(create_job, setups; root=joinpath(tmpdir, "study"))

        @test length(study.jobs) == 1
        @test study.jobs[1].spatial_setup isa MultibodySetup

        submit!(study)
        @test study.sim_success[1] == true
    end

    Peridynamics.MPI_RUN[] = mpi_run_current_value
end

@testitem "Study logfile format and content" tags=[:skipci] begin
    mpi_run_current_value = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

    mktempdir() do tmpdir
        function create_job(setup::NamedTuple, root::String)
            body = Body(BBMaterial(), rand(3, 10), rand(10))
            material!(body, horizon=1, E=setup.E, rho=1, Gc=1)
            velocity_ic!(body, :all_points, :x, 1.0)
            vv = VelocityVerlet(steps=setup.n_steps)
            path = joinpath(root, "sim_E$(Int(setup.E))")
            job = Job(body, vv; path=path, freq=5)
            return job
        end

        setups = [
            (; E=1e9, n_steps=5),
            (; E=2e9, n_steps=10),
        ]

        study = Study(create_job, setups; root=joinpath(tmpdir, "study"))
        submit!(study; quiet=true)

        @test isfile(study.logfile)
        logcontent = read(study.logfile, String)

        # Check header
        @test contains(logcontent, "SIMULATION STUDY LOGFILE")

        # Check all parameters are logged
        @test contains(logcontent, "E: 1.0e9")
        @test contains(logcontent, "E: 2.0e9")
        @test contains(logcontent, "n_steps: 5")
        @test contains(logcontent, "n_steps: 10")

        # Check paths
        @test contains(logcontent, "sim_E1000000000")
        @test contains(logcontent, "sim_E2000000000")

        # Check status
        @test count("completed ✓", logcontent) == 2

        # Check timing info present
        @test contains(logcontent, "seconds")
    end

    Peridynamics.MPI_RUN[] = mpi_run_current_value
end
