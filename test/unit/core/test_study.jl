# `Study`: a set of jobs created from setups, submitted in one go and post-processed per job
# (`src/core/study.jl`). All items build their jobs with `Fixtures.tiny_job` (a ten-point line
# body, 5–15 explicit steps) so the study machinery is compiled once for the whole file.

@testitem "Study: constructor" setup=[Fixtures] begin
    root = joinpath(mktempdir(), "study")
    study = Study(Fixtures.tiny_job, Fixtures.TINY_SETUPS; root)
    @test study isa Study
    @test length(study.jobs) == length(study.setups) == length(study.jobpaths) == 3
    @test all(.!study.sim_success)
    @test study.root == root
    @test study.logfile == joinpath(root, "study_log.log")
    @test all(contains.(study.jobpaths, "sim_"))
    # a custom logfile name
    study2 = Study(Fixtures.tiny_job, Fixtures.TINY_SETUPS; root, logfile_name="my_log.log")
    @test study2.logfile == joinpath(root, "my_log.log")
end

@testitem "Study: constructor errors" setup=[Fixtures] begin
    root = mktempdir()
    # no setups
    @test_throws ArgumentError Study(Fixtures.tiny_job, NamedTuple[]; root)
    # setups with different fields
    setups = [(; E=1.0, n_steps=5, velocity=1.0), (; E=1.0, n_steps=5, force=2.0)]
    @test_throws ArgumentError Study(Fixtures.tiny_job, setups; root)
    # a job creator that errors for a setup (negative number of steps)
    setups = [(; E=1.0, n_steps=5), (; E=1.0, n_steps=-1)]
    @test_throws ArgumentError Study(Fixtures.tiny_job, setups; root)
    # two setups that map to the same job path
    function same_path_job(setup, root)
        body = Fixtures.line10(BBMaterial(); Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        return Job(body, VelocityVerlet(steps=setup.n_steps); path=joinpath(root, "sim"))
    end
    setups = [(; n_steps=10), (; n_steps=20)]
    @test_throws ArgumentError Study(same_path_job, setups; root)
end

@testitem "Study: show" setup=[Fixtures] begin
    study = Study(Fixtures.tiny_job, Fixtures.TINY_SETUPS[1:2]; root=mktempdir())
    msg = sprint(show, MIME("text/plain"), study)
    @test contains(msg, "Study with 2 jobs")
    @test contains(msg, "✗") # not yet submitted
    @test contains(msg, "sim_1.0_5")
    @test contains(msg, "sim_2.0_10")
end

@testitem "submit!: all jobs succeed" tags=[:simulation] setup=[Fixtures] begin
    root = joinpath(mktempdir(), "my_study")
    study = Study(Fixtures.tiny_job, Fixtures.TINY_SETUPS[1:2]; root)
    @test !isdir(root)
    submit!(study; quiet=true)
    @test all(study.sim_success)
    # the directory structure: the study root, the logfile, one directory per job
    @test isdir(root)
    @test isfile(study.logfile)
    @test all(isdir, study.jobpaths)
    @test all(isfile(joinpath(path, "logfile.log")) for path in study.jobpaths)
    # the logfile: header, parameters, paths, status, timing, numbering
    logcontent = read(study.logfile, String)
    @test contains(logcontent, "SIMULATION STUDY LOGFILE")
    @test contains(logcontent, "E: 1.0") && contains(logcontent, "E: 2.0")
    @test contains(logcontent, "n_steps: 5") && contains(logcontent, "n_steps: 10")
    @test contains(logcontent, "sim_1.0_5") && contains(logcontent, "sim_2.0_10")
    @test count("completed ✓", logcontent) == 2
    @test contains(logcontent, "seconds")
    @test contains(logcontent, "(1/2) Simulation") && contains(logcontent, "(2/2) Simulation")
end

@testitem "submit!: failing jobs are recorded and the others still run" tags=[:simulation] setup=[Fixtures] begin
    setups = [(; E=1.0, n_steps=5, fail=false), (; E=2.0, n_steps=10, fail=true),
              (; E=3.0, n_steps=15, fail=false)]
    study = Study(Fixtures.tiny_job, setups; root=joinpath(mktempdir(), "study"))
    submit!(study; quiet=true)
    @test study.sim_success == [true, false, true]
    logcontent = read(study.logfile, String)
    @test count("completed ✓", logcontent) == 2
    @test count("failed ✗", logcontent) == 1
    @test contains(logcontent, "(1/3) Simulation") && contains(logcontent, "(3/3) Simulation")
    # the error of the failing job is in its own logfile
    @test contains(read(joinpath(study.jobpaths[2], "logfile.log"), String),
                   "deliberately failing job")
end

@testitem "submit!: the loud path prints and a failing job does not crash it" tags=[:simulation] setup=[Fixtures] begin
    # `quiet=false` is the default of `submit!`; the study progress and the error of the
    # failing job are printed, which is captured here and must not throw
    setups = [(; E=1.0, n_steps=5, fail=false), (; E=2.0, n_steps=5, fail=true)]
    study = Study(Fixtures.tiny_job, setups; root=joinpath(mktempdir(), "study"))
    printed = mktemp() do path, io
        redirect_stdio(; stdout=io, stderr=io) do
            submit!(study)
        end
        flush(io)
        read(path, String)
    end
    Peridynamics.set_quiet!(true) # `submit` left the global quiet flag cleared
    @test contains(printed, "MULTITHREADING SIMULATION") # the banner of each job
    @test contains(printed, "failed with error")
    @test contains(printed, "Study completed (1 / 2 simulations successful)")
    @test study.sim_success == [true, false]
end

@testitem "submit!: resuming skips the jobs a previous run completed" tags=[:simulation] setup=[Fixtures] begin
    setups = [(; E=1.0, n_steps=1), (; E=2.0, n_steps=1)]
    study = Study(Fixtures.tiny_job, setups; root=joinpath(mktempdir(), "study"))
    # an interrupted run: the study root exists and the logfile marks the first job completed
    mkpath(study.root)
    mkpath(study.jobpaths[1])
    open(study.logfile, "w") do io
        write(io, "SIMULATION STUDY LOGFILE\n\n")
        write(io, "Simulation `$(study.jobpaths[1])`:\n")
        write(io, "  status: completed ✓ (0.00 seconds)\n\n")
    end
    submit!(study; quiet=true)
    @test study.sim_success == [true, true]
    @test isdir(study.jobpaths[2])
    logcontent = read(study.logfile, String)
    @test contains(logcontent, "--- RESUMED:")
    @test contains(logcontent, "status: skipped (completed in a previous run)")
    # the processing recognizes the completed job as well
    processed = process_each_job((job, setup) -> (; E=setup.E), study, (; E=0.0))
    @test processed[1].E == 1.0 && processed[2].E == 2.0
end

@testitem "Study: with a MultibodySetup" tags=[:simulation] setup=[Fixtures] begin
    function multibody_job(setup, root)
        b1 = Fixtures.line10(BBMaterial(); Gc=1)
        velocity_ic!(b1, :all_points, :x, setup.velocity)
        b2 = Fixtures.line10(OSBMaterial(); Gc=1)
        b2.position[2, :] .+= 2.0 # next to the first body
        ms = MultibodySetup(:body1 => b1, :body2 => b2)
        return Job(ms, VelocityVerlet(steps=setup.n_steps); path=joinpath(root, "sim_$(setup.id)"), freq=5)
    end
    study = Study(multibody_job, [(; id=1, n_steps=5, velocity=1.0)]; root=joinpath(mktempdir(), "study"))
    @test study.jobs[1].spatial_setup isa MultibodySetup
    submit!(study; quiet=true)
    @test study.sim_success[1]
end

@testitem "process_each_job: results, failed jobs and processing errors" tags=[:simulation] setup=[Fixtures] begin
    setups = [(; E=1.0, n_steps=5, fail=false), (; E=2.0, n_steps=10, fail=true),
              (; E=3.0, n_steps=15, fail=false)]
    study = Study(Fixtures.tiny_job, setups; root=joinpath(mktempdir(), "study"))
    submit!(study; quiet=true)
    @test study.sim_success == [true, false, true]

    # the processing function sees the job and the setup; failed jobs get the default result
    default_result = (; E=0.0, n_steps=0, n_files=0)
    process(job, setup) = (; setup.E, setup.n_steps, n_files=length(readdir(job.options.root)))
    results = @test_logs (:warn, r"skipping processing for failed job") process_each_job(process, study, default_result)
    @test length(results) == 3
    @test results[1] == (; E=1.0, n_steps=5, n_files=results[1].n_files) && results[1].n_files > 0
    @test results[2] == default_result
    @test results[3].E == 3.0 && results[3].n_steps == 15 && results[3].n_files > 0

    # `process_failed=true` processes the failed job as well
    results = process_each_job(process, study, default_result; process_failed=true)
    @test results[2].E == 2.0 && results[2].n_steps == 10

    # without MPI `only_root` changes nothing
    counter = Ref(0)
    counting(job, setup) = (counter[] += 1; (; setup.E, n_steps=0, n_files=0))
    results = @test_logs (:warn, r"skipping processing") process_each_job(counting, study, default_result; only_root=true)
    @test counter[] == 2
    results = process_each_job(counting, study, default_result; only_root=true, process_failed=true)
    @test counter[] == 5

    # a processing error gives the default result, is printed to stderr and written to an
    # error logfile in the job directory
    erroring(job, setup) = setup.E == 3.0 ? error("intentional processing error") : (; setup.E, n_steps=0, n_files=0)
    results, printed = mktemp() do path, io
        r = redirect_stderr(io) do
            @test_logs (:warn, r"skipping processing") process_each_job(erroring, study, default_result)
        end
        flush(io)
        r, read(path, String)
    end
    @test contains(printed, "ERROR:") && contains(printed, "intentional processing error")
    @test results[1].E == 1.0
    @test results[3] == default_result
    error_logs = filter(f -> occursin("proc_error.log", f), readdir(study.jobpaths[3]))
    @test length(error_logs) == 1
    error_log = read(joinpath(study.jobpaths[3], error_logs[1]), String)
    @test contains(error_log, "ERROR:") && contains(error_log, "intentional processing error")

    # an empty result type is fine
    empty = @test_logs (:warn, r"skipping processing") process_each_job((job, setup) -> NamedTuple(), study, NamedTuple())
    @test empty == [NamedTuple(), NamedTuple(), NamedTuple()]
end
