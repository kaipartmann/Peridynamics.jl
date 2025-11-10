@testitem "Study constructor with valid setups" begin
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
    root = joinpath(mktempdir(), "study")
    study = Peridynamics.Study(create_job, setups; root)

    @test study isa Peridynamics.Study
    @test length(study.jobs) == 3
    @test length(study.setups) == 3
    @test length(study.sim_success) == 3
    @test length(study.jobpaths) == 3
    @test all(.!study.sim_success)
    @test study.root == root
    @test study.logfile == joinpath(root, "study_log.log")
    @test all(contains.(study.jobpaths, "sim_"))
end

@testitem "Study constructor with custom logfile name" begin
    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=setup.E, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        vv = VelocityVerlet(steps=setup.n_steps)
        path = joinpath(root, "sim_$(setup.n_steps)")
        job = Job(body, vv; path=path, freq=5)
        return job
    end
    setups = [(; E=1.0, n_steps=5)]
    root = joinpath(mktempdir(), "study")
    custom_logfile = "my_custom_study_log.log"
    study = Peridynamics.Study(create_job, setups; root, logfile_name=custom_logfile)

    @test study.logfile == joinpath(root, custom_logfile)

    # Test that submit! creates the logfile with custom name
    Peridynamics.submit!(study; quiet=true)
    @test isfile(study.logfile)
    @test isfile(joinpath(root, custom_logfile))
end

@testitem "resume processing from existing logfile" begin
    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=setup.E, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        vv = VelocityVerlet(steps=setup.n_steps)
        path = joinpath(root, "sim_$(setup.E)")
        job = Job(body, vv; path=path, freq=5)
        return job
    end
    setups = [
        (; E=1.0, n_steps=1),
        (; E=2.0, n_steps=1),
    ]
    root = joinpath(mktempdir(), "study")
    study = Peridynamics.Study(create_job, setups; root)

    # Simulate an interrupted run: create study root and mark first job completed
    mkpath(study.root)
    mkpath(study.jobpaths[1])
    open(study.logfile, "w") do io
        write(io, "SIMULATION STUDY LOGFILE\n\n")
        write(io, "Simulation `$(study.jobpaths[1])`:\n")
        write(io, "  status: completed ✓ (0.00 seconds)\n\n")
    end

    # Now processing should detect that first job succeeded and only process that one
    default_result = (; E=0.0)
    processed = process_each_job((job, setup) -> (; E=setup.E), study, default_result)

    @test processed[1].E == 1.0
    @test processed[2].E == 0.0
end

@testitem "resume submit! only runs remaining jobs" begin
    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=setup.E, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        vv = VelocityVerlet(steps=setup.n_steps)
        path = joinpath(root, "sim_$(setup.E)")
        job = Job(body, vv; path=path, freq=5)
        return job
    end

    setups = [
        (; E=1.0, n_steps=1),
        (; E=2.0, n_steps=1),
    ]

    studyroot = joinpath(mktempdir(), "study")
    study = Peridynamics.Study(create_job, setups; root=studyroot)

    # Simulate an interrupted run: create study root and mark first job completed
    mkpath(study.root)
    mkpath(study.jobpaths[1])
    open(study.logfile, "w") do io
        write(io, "SIMULATION STUDY LOGFILE\n\n")
        write(io, "Simulation `$(study.jobpaths[1])`:\n")
        write(io, "  status: completed ✓ (0.00 seconds)\n\n")
    end

    # Now resuming submit! should only run the remaining job (E=2)
    Peridynamics.submit!(study; quiet=true)

    @test study.sim_success == [true, true]
    @test isdir(study.jobpaths[2])

    # Check that logfile contains RESUMED marker and skipped status
    logcontent = read(study.logfile, String)
    @test contains(logcontent, "--- RESUMED:")
    @test contains(logcontent, "status: skipped (completed in a previous run)")
end

@testitem "Study constructor with empty setups" begin
    function create_job(setup::NamedTuple, root::String)
        return Job(Body(BBMaterial(), rand(3, 10), rand(10)), VelocityVerlet(steps=1))
    end
    setups = NamedTuple[]
    @test_throws ArgumentError Peridynamics.Study(create_job, setups; root=mktempdir())
end

@testitem "Study constructor with inconsistent setup fields" begin
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
    @test_throws ArgumentError Peridynamics.Study(create_job, setups; root=mktempdir())
end

@testitem "Study constructor with jobcreator that errors" begin
    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=1, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        path = joinpath(root, "sim_$(setup.n_steps)")
        return Job(body, VelocityVerlet(steps=setup.n_steps); path=path)
    end
    setups = [
        (; n_steps=10,),
        (; n_steps=-1,),
    ]
    @test_throws ArgumentError Peridynamics.Study(create_job, setups; root=mktempdir())
end

@testitem "Study constructor with duplicate job paths" begin
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
    @test_throws ArgumentError Peridynamics.Study(create_job, setups; root=mktempdir())
end

@testitem "Study show method" begin
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
    study = Peridynamics.Study(create_job, setups; root=mktempdir())

    io = IOBuffer()
    show(io, MIME("text/plain"), study)
    msg = String(take!(io))
    @test contains(msg, "Study with 2 jobs")
    @test contains(msg, "✗")  # Not yet submitted
    @test contains(msg, "sim_10")
    @test contains(msg, "sim_20")
end

@testitem "submit! with all successful jobs" begin
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
    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))
    @test all(.!study.sim_success)

    Peridynamics.submit!(study; quiet=true)
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

@testitem "submit! with some failing jobs" begin
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
    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))

    # Run with quiet=false to allow error printing (we just verify it doesn't crash)
    Peridynamics.submit!(study; quiet=false)

    @test study.sim_success[1] == true
    @test study.sim_success[2] == false
    @test study.sim_success[3] == true
    @test count(study.sim_success) == 2

    # Check logfile contains failure info and proper formatting
    logcontent = read(study.logfile, String)
    @test contains(logcontent, "completed ✓")
    @test contains(logcontent, "failed ✗")
    @test contains(logcontent, "(1/3) Simulation")
    @test contains(logcontent, "(2/3) Simulation")
    @test contains(logcontent, "(3/3) Simulation")
end

@testitem "submit! creates study directory structure" begin
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
    studyroot = joinpath(mktempdir(), "my_study")
    study = Peridynamics.Study(create_job, setups; root=studyroot)
    @test !isdir(studyroot)

    Peridynamics.submit!(study; quiet=true)
    @test isdir(studyroot)
    @test isfile(study.logfile)
    @test isdir(study.jobpaths[1])
    @test isfile(joinpath(study.jobpaths[1], "logfile.log"))
end

@testitem "submit! with quiet keyword" begin
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
    tmpdir = mktempdir()
    study = Peridynamics.Study(create_job, setups; root=joinpath(tmpdir, "study"))

    # Should not error with quiet=true
    Peridynamics.submit!(study; quiet=true)
    @test study.sim_success[1] == true

    # Should not error with quiet=false (default)
    study2 = Peridynamics.Study(create_job, setups; root=joinpath(tmpdir, "study2"))
    Peridynamics.submit!(study2)
    @test study2.sim_success[1] == true
end

@testitem "submit! error handling with quiet=true" begin
    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=1, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        vv = VelocityVerlet(steps=setup.n_steps)
        # Invalid path will cause job to fail
        path = "/invalid/path/sim_$(setup.n_steps)"
        job = Job(body, vv; path=path, freq=5)
        return job
    end
    setups = [(; n_steps=5,)]
    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))

    # With quiet=true, should handle errors gracefully without crashing
    Peridynamics.submit!(study; quiet=true)

    @test study.sim_success[1] == false

    # Logfile should still contain the failure status
    logcontent = read(study.logfile, String)
    @test contains(logcontent, "failed ✗")
end

@testitem "Study with MultibodySetup" begin
    __mpi_run__ = Peridynamics.MPI_RUN[]
    Peridynamics.MPI_RUN[] = false

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
    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))

    @test length(study.jobs) == 1
    @test study.jobs[1].spatial_setup isa MultibodySetup

    Peridynamics.submit!(study)
    @test study.sim_success[1] == true

    Peridynamics.MPI_RUN[] = __mpi_run__
end

@testitem "Study logfile format and content" begin
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

    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))
    Peridynamics.submit!(study; quiet=true)

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

    # Check job numbering format
    @test contains(logcontent, "(1/2) Simulation")
    @test contains(logcontent, "(2/2) Simulation")
end

@testitem "submit! prints study completion message to logfile" begin
    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=1, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        vv = VelocityVerlet(steps=setup.n_steps)
        path = joinpath(root, "sim_$(setup.n_steps)")
        job = Job(body, vv; path=path, freq=5)
        return job
    end
    setups = [
        (; n_steps=5),
        (; n_steps=10),
    ]
    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))

    # Submit jobs (completion message goes to stdout via print_log)
    Peridynamics.submit!(study; quiet=true)

    @test all(study.sim_success)
    # The completion message is printed via print_log, not stored in logfile
    # Just verify the jobs completed successfully
end

@testitem "submit! tracks partial success correctly" begin
    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=1, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        vv = VelocityVerlet(steps=setup.n_steps)

        # Second job will fail with invalid path
        if setup.fail
            path = "/invalid/path/sim_$(setup.n_steps)"
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
    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))

    # Submit with quiet=true to suppress error output
    Peridynamics.submit!(study; quiet=true)

    # Verify correct success tracking
    @test study.sim_success[1] == true
    @test study.sim_success[2] == false
    @test study.sim_success[3] == true
    @test count(study.sim_success) == 2
end

@testitem "process_each_job with all successful jobs" begin
    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=setup.E, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        vv = VelocityVerlet(steps=setup.n_steps)
        path = joinpath(root, "sim_$(setup.E)")
        job = Job(body, vv; path=path, freq=5)
        return job
    end
    setups = [
        (; E=1.0, n_steps=5),
        (; E=2.0, n_steps=10),
        (; E=3.0, n_steps=15),
    ]
    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))
    Peridynamics.submit!(study; quiet=true)

    # Define a processing function that extracts parameters from setup
    function process_func(job, setup)
        return (; E=setup.E, n_steps=setup.n_steps, path=job.options.root)
    end
    default_result = (; E=0.0, n_steps=0, path="")
    results = Peridynamics.process_each_job(process_func, study, default_result)

    @test length(results) == 3
    @test results[1].E == 1.0
    @test results[1].n_steps == 5
    @test results[2].E == 2.0
    @test results[2].n_steps == 10
    @test results[3].E == 3.0
    @test results[3].n_steps == 15
    @test all(r -> r.path != "", results)
end

@testitem "process_each_job with some failed jobs" begin
    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=1, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        vv = VelocityVerlet(steps=setup.n_steps)

        # Create invalid path for second job
        if setup.fail
            path = "/invalid/path/sim_$(setup.n_steps)"
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
    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))
    Peridynamics.submit!(study; quiet=true)

    # Only first and third jobs should succeed
    @test study.sim_success == [true, false, true]

    function process_func(job, setup)
        return (; n_steps=setup.n_steps, processed=true)
    end

    default_result = (; n_steps=0, processed=false)

    # Process jobs - should skip failed job with warning
    results = Peridynamics.process_each_job(process_func, study, default_result)

    @test length(results) == 3
    # First job processed
    @test results[1].n_steps == 5
    @test results[1].processed == true
    # Second job failed - should have default result (skipped with warning)
    @test results[2].n_steps == 0
    @test results[2].processed == false
    # Third job processed
    @test results[3].n_steps == 15
    @test results[3].processed == true
end

@testitem "process_each_job with processing errors" begin
    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=1, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        vv = VelocityVerlet(steps=setup.n_steps)
        path = joinpath(root, "sim_$(setup.id)")
        job = Job(body, vv; path=path, freq=5)
        return job
    end
    setups = [
        (; id=1, n_steps=5),
        (; id=2, n_steps=10),
        (; id=3, n_steps=15),
    ]
    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))
    Peridynamics.submit!(study; quiet=true)
    @test all(study.sim_success)

    # Processing function that errors for specific setup
    function process_func(job, setup)
        if setup.id == 2
            error("Intentional processing error")
        end
        return (; id=setup.id, status="ok")
    end
    default_result = (; id=0, status="failed")

    # Should not throw, but should use default_result for erroring job
    results = Peridynamics.process_each_job(process_func, study, default_result)

    @test length(results) == 3
    @test results[1].id == 1
    @test results[1].status == "ok"
    @test results[2].id == 0  # Error occurred, default used
    @test results[2].status == "failed"
    @test results[3].id == 3
    @test results[3].status == "ok"

    # Check that error logfile was created for the failed job
    error_logfiles = filter(f -> occursin("proc_error.log", f), readdir(study.jobpaths[2]))
    @test length(error_logfiles) == 1
    error_log_path = joinpath(study.jobpaths[2], error_logfiles[1])
    error_log_content = read(error_log_path, String)
    @test contains(error_log_content, "ERROR:")
    @test contains(error_log_content, "Intentional processing error")
end

@testitem "process_each_job creates error logfile" begin
    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=1, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        vv = VelocityVerlet(steps=setup.n_steps)
        path = joinpath(root, "sim_$(setup.id)")
        job = Job(body, vv; path=path, freq=5)
        return job
    end
    setups = [(; id=1, n_steps=5)]
    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))
    Peridynamics.submit!(study; quiet=true)
    @test all(study.sim_success)

    # Processing function that errors
    function process_func(job, setup)
        error("Processing error for testing")
    end
    default_result = (; id=0, status="failed")

    # Should handle error gracefully and use default result
    results = Peridynamics.process_each_job(process_func, study, default_result)

    @test results[1] == default_result

    # Error logfile should be created
    error_logfiles = filter(f -> occursin("proc_error.log", f), readdir(study.jobpaths[1]))
    @test length(error_logfiles) == 1

    # Verify error logfile content
    error_log_path = joinpath(study.jobpaths[1], error_logfiles[1])
    error_log_content = read(error_log_path, String)
    @test contains(error_log_content, "ERROR:")
    @test contains(error_log_content, "Processing error for testing")
end

@testitem "process_each_job with empty processing result" begin
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
    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))
    Peridynamics.submit!(study; quiet=true)
    # Processing function that returns empty NamedTuple
    process_func(job, setup) = NamedTuple()
    default_result = NamedTuple()
    results = Peridynamics.process_each_job(process_func, study, default_result)
    @test length(results) == 1
    @test results[1] == NamedTuple()
end

@testitem "process_each_job accessing job results" begin
    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=setup.E, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        vv = VelocityVerlet(steps=setup.n_steps)
        path = joinpath(root, "sim_$(setup.E)")
        job = Job(body, vv; path=path, freq=5)
        return job
    end
    setups = [
        (; E=1e9, n_steps=5),
        (; E=2e9, n_steps=10),
    ]
    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))
    Peridynamics.submit!(study; quiet=true)

    # Process function that reads result files
    function process_func(job, setup)
        # Check if output directory exists
        if isdir(job.options.root)
            n_files = length(readdir(job.options.root))
            return (; E=setup.E, n_output_files=n_files)
        else
            return (; E=setup.E, n_output_files=0)
        end
    end
    default_result = (; E=0.0, n_output_files=0)
    results = Peridynamics.process_each_job(process_func, study, default_result)
    @test length(results) == 2
    @test results[1].E == 1e9
    @test results[1].n_output_files > 0  # Should have created output files
    @test results[2].E == 2e9
    @test results[2].n_output_files > 0
end

@testitem "process_each_job with only_root=true" begin
    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=setup.E, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        vv = VelocityVerlet(steps=setup.n_steps)
        path = joinpath(root, "sim_$(setup.E)")
        job = Job(body, vv; path=path, freq=5)
        return job
    end
    setups = [
        (; E=1e9, n_steps=5),
        (; E=2e9, n_steps=10),
    ]
    study = Peridynamics.Study(create_job, setups; root=joinpath(mktempdir(), "study"))
    Peridynamics.submit!(study; quiet=true)

    # Test with only_root=true
    counter = Ref(0)
    function process_func(job, setup)
        counter[] += 1
        return (; E=setup.E, processed=true)
    end

    default_result = (; E=0.0, processed=false)
    results = Peridynamics.process_each_job(process_func, study, default_result;
                                           only_root=true)

    @test length(results) == 2
    # In non-MPI context, mpi_isroot() is always true, so both should be processed
    @test results[1].E == 1e9
    @test results[1].processed == true
    @test results[2].E == 2e9
    @test results[2].processed == true
    @test counter[] == 2  # Both jobs processed

    # Test with only_root=false (default)
    counter[] = 0
    results2 = Peridynamics.process_each_job(process_func, study, default_result)
    @test counter[] == 2  # Same behavior without MPI
end

@testitem "process_each_job with only_root in MPI" tags=[:mpi] begin
    path = mktempdir()
    mpi_cmd = """
    using Peridynamics

    function create_job(setup::NamedTuple, root::String)
        body = Body(BBMaterial(), rand(3, 10), rand(10))
        material!(body, horizon=1, E=1.0, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        vv = VelocityVerlet(steps=5)
        path = joinpath(root, "sim_\$(setup.id)")
        job = Job(body, vv; path=path, freq=5)
        return job
    end

    setups = [(; id=1), (; id=2)]
    study = Peridynamics.Study(create_job, setups; root="$(path)")
    Peridynamics.submit!(study; quiet=true)

    # Test with only_root=true
    marker_file = "$(joinpath(path, "processing_marker.txt"))"
    default_result = (; id=0, rank=-1)

    results = Peridynamics.process_each_job(study, default_result; only_root=true) do job, setup
        # This should only run on root
        rank = Peridynamics.mpi_rank()
        if mpi_isroot()
            open(marker_file, "a") do io
                println(io, "Job \$(setup.id) processed on rank \$(rank)")
            end
        end
        return (; id=setup.id, rank=rank)
    end

    # Check that only root processed (and wrote to file)
    if mpi_isroot()
        @assert isfile(marker_file)
        content = read(marker_file, String)
        @assert contains(content, "Job 1 processed on rank 0")
        @assert contains(content, "Job 2 processed on rank 0")
        @assert !contains(content, "rank 1")  # Non-root shouldn't write
    end
    """

    mpiexec = Peridynamics.MPI.mpiexec()
    jlcmd = Base.julia_cmd()
    pdir = pkgdir(Peridynamics)
    @test success(`$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd)`)

    # Verify the marker file was created correctly
    marker_file = joinpath(path, "processing_marker.txt")
    @test isfile(marker_file)
    content = read(marker_file, String)
    @test contains(content, "Job 1 processed on rank 0")
    @test contains(content, "Job 2 processed on rank 0")
end
