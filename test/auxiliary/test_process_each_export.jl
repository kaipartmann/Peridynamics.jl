@testitem "process_each_export - basic functionality" begin
    root = mktempdir()
    root_post_threads = joinpath(root, "post_threads")
    mkpath(root_post_threads)
    root_post_serial = joinpath(root, "post_serial")
    mkpath(root_post_serial)

    l, Δx, δ, a = 1.0, 1 / 4, 3.015 / 4, 0.5
    pos, vol = uniform_box(l, l, l, Δx)
    b1 = Body(BBMaterial(), pos, vol)
    material!(b1; horizon=3.015Δx, E=2.1e5, rho=8e-6)
    point_set!(y -> y > l / 2 - Δx, b1, :set_top)
    point_set!(y -> y < -l / 2 + Δx, b1, :set_bottom)
    velocity_bc!(t -> 30, b1, :set_top, :y)
    velocity_bc!(t -> -30, b1, :set_bottom, :y)
    vv = VelocityVerlet(steps=2)
    job = Job(b1, vv; path=root, freq=1)

    @test_throws ArgumentError process_each_export((r0, r, id) -> nothing, job)

    submit(job)

    # Test threads backend (legacy mode without result collection)
    process_each_export(job) do result0, result, file_id
        filename = string("max_displacement_", file_id, ".txt")
        open(joinpath(root_post_threads, filename), "w+") do io
            maxdisp = maximum(result[:displacement][1,:])
            msg = string("maximum displacement x: ", maxdisp)
            write(io, msg)
        end
        return nothing
    end
    file_1_threads = joinpath(root_post_threads, "max_displacement_1.txt")
    @test isfile(file_1_threads)
    @test contains(read(file_1_threads, String), "maximum displacement x: 0.0")
    file_2_threads = joinpath(root_post_threads, "max_displacement_2.txt")
    @test isfile(file_2_threads)
    @test contains(read(file_2_threads, String), "maximum displacement x: 0.0")
    file_3_threads = joinpath(root_post_threads, "max_displacement_3.txt")
    @test isfile(file_3_threads)
    @test contains(read(file_3_threads, String), "maximum displacement x: 2.4")

    # Test serial backend (legacy mode without result collection)
    process_each_export(job; serial=true) do result0, result, file_id
        filename = string("max_displacement_", file_id, ".txt")
        open(joinpath(root_post_serial, filename), "w+") do io
            msg = string("maximum displacement x: ", maximum(result[:displacement][1,:]))
            write(io, msg)
        end
        return nothing
    end
    file_1_serial = joinpath(root_post_serial, "max_displacement_1.txt")
    @test isfile(file_1_serial)
    @test contains(read(file_1_serial, String), "maximum displacement x: 0.0")
    file_2_serial = joinpath(root_post_serial, "max_displacement_2.txt")
    @test isfile(file_2_serial)
    @test contains(read(file_2_serial, String), "maximum displacement x: 0.0")
    file_3_serial = joinpath(root_post_serial, "max_displacement_3.txt")
    @test isfile(file_3_serial)
    @test contains(read(file_3_serial, String), "maximum displacement x: 2.4")
end

@testitem "find_vtk_files" begin
    root = mktempdir()
    vtk_files_unsorted = ["timestep_123.pvtu", "abcd_timestep_000005.pvtu",
                          "timestep_02.pvtu", "_1.pvtu"]
    for file in vtk_files_unsorted
        open(joinpath(root, file), "w+") do io
            write(io, "")
        end
    end

    vtk_files = Peridynamics.find_vtk_files(root)
    @test basename.(vtk_files) == ["_1.pvtu", "timestep_02.pvtu",
                                   "abcd_timestep_000005.pvtu", "timestep_123.pvtu"]
end

@testitem "process_each_export with barrier" begin
    root = mktempdir()
    l, Δx = 1.0, 1 / 4
    pos, vol = uniform_box(l, l, l, Δx)
    b1 = Body(BBMaterial(), pos, vol)
    material!(b1; horizon=3.015Δx, E=2.1e5, rho=8e-6)
    point_set!(y -> y > l / 2 - Δx, b1, :set_top)
    point_set!(y -> y < -l / 2 + Δx, b1, :set_bottom)
    velocity_bc!(t -> 30, b1, :set_top, :y)
    velocity_bc!(t -> -30, b1, :set_bottom, :y)
    vv = VelocityVerlet(steps=2)
    job = Job(b1, vv; path=root, freq=1)
    submit(job)

    # Test that barrier parameter works without error in non-MPI context
    counter = Ref(0)
    process_each_export(job; serial=true, barrier=false) do r0, r, id
        counter[] += 1
    end
    @test counter[] == 3  # reference + 2 time steps

    # Test with barrier=true (should also work in non-MPI context)
    counter[] = 0
    process_each_export(job; serial=true, barrier=true) do r0, r, id
        counter[] += 1
    end
    @test counter[] == 3
end

@testitem "process_each_export with result collection - threads" begin
    # initial setup
    Peridynamics.MPI_RUN[] = false

    using Base.Threads: nthreads

    root = mktempdir()
    l, Δx = 1.0, 1 / 4
    pos, vol = uniform_box(l, l, l, Δx)
    b1 = Body(BBMaterial(), pos, vol)
    material!(b1; horizon=3.015Δx, E=2.1e5, rho=8e-6)
    point_set!(y -> y > l / 2 - Δx, b1, :set_top)
    point_set!(y -> y < -l / 2 + Δx, b1, :set_bottom)
    velocity_bc!(t -> 30, b1, :set_top, :y)
    velocity_bc!(t -> -30, b1, :set_bottom, :y)
    vv = VelocityVerlet(steps=5)  # More steps for better testing
    job = Job(b1, vv; path=root, freq=1)
    submit(job)

    # Test result collection with NamedTuples (use threads if serial=false and !mpi_run())
    default_value = (; max_disp=NaN, file_id=0)
    results = process_each_export(job, default_value; serial=false) do r0, r, id
        max_disp = maximum(r[:displacement])
        return (; max_disp, file_id=id)
    end

    @test results isa Vector{NamedTuple{(:max_disp, :file_id), Tuple{Float64, Int}}}
    @test length(results) == 6  # reference + 5 time steps
    # Verify all file IDs are present and correct
    @test all(results[i].file_id == i for i in 1:6)
    @test results[1].max_disp ≈ 0.0
    @test results[end].max_disp > 0.0  # Final step should have displacement

    # Test without result collection (legacy mode) - just verify it returns nothing
    result_nothing = process_each_export(job; serial=false) do r0, r, id
        return nothing
    end
    @test result_nothing === nothing

    # Test that results are computed correctly in parallel, each thread calcs independently
    results_parallel = process_each_export(job, 0.0; serial=false) do r0, r, id
        return id * 2.5
    end
    @test results_parallel == [2.5, 5.0, 7.5, 10.0, 12.5, 15.0]

    # Test that NamedTuples also work correctly
    default_value = (; id=0, squared=0)
    results_parallel = process_each_export(job, default_value; serial=false) do r0, r, id
        return (; id=id, squared=id*id)
    end
    @test all(results_parallel[i].id == i for i in 1:6)
    @test all(results_parallel[i].squared == i*i for i in 1:6)
end

@testitem "process_each_export with result collection - serial" begin
    root = mktempdir()
    l, Δx = 1.0, 1 / 4
    pos, vol = uniform_box(l, l, l, Δx)
    b1 = Body(BBMaterial(), pos, vol)
    material!(b1; horizon=3.015Δx, E=2.1e5, rho=8e-6)
    point_set!(y -> y > l / 2 - Δx, b1, :set_top)
    point_set!(y -> y < -l / 2 + Δx, b1, :set_bottom)
    velocity_bc!(t -> 30, b1, :set_top, :y)
    velocity_bc!(t -> -30, b1, :set_bottom, :y)
    vv = VelocityVerlet(steps=2)
    job = Job(b1, vv; path=root, freq=1)
    submit(job)

    # Test result collection in serial mode
    default_value = (; avg_ux=NaN, avg_uy=NaN)
    results = process_each_export(job, default_value; serial=true) do r0, r, id
        ux = @view r[:displacement][1, :]
        uy = @view r[:displacement][2, :]
        avg_ux = sum(ux) / length(ux)
        avg_uy = sum(uy) / length(uy)
        return (; avg_ux, avg_uy)
    end

    @test results isa Vector{NamedTuple{(:avg_ux, :avg_uy), Tuple{Float64, Float64}}}
    @test length(results) == 3
    @test all(.!isnan.(getfield.(results, :avg_ux)))
    @test all(.!isnan.(getfield.(results, :avg_uy)))
    @test results[1].avg_ux ≈ 0.0
    @test results[1].avg_uy ≈ 0.0
end

@testitem "process_each_export - MPI tests" tags=[:mpi] begin
    # Run simulation once and use for all MPI tests
    root = mktempdir()
    l, Δx = 1.0, 1 / 4
    pos, vol = uniform_box(l, l, l, Δx)
    b1 = Body(BBMaterial(), pos, vol)
    material!(b1; horizon=3.015Δx, E=2.1e5, rho=8e-6)
    point_set!(y -> y > l / 2 - Δx, b1, :set_top)
    point_set!(y -> y < -l / 2 + Δx, b1, :set_bottom)
    velocity_bc!(t -> 30, b1, :set_top, :y)
    velocity_bc!(t -> -30, b1, :set_bottom, :y)
    vv = VelocityVerlet(steps=5)
    job = Job(b1, vv; path=root, freq=1)
    submit(job)

    mpiexec = Peridynamics.MPI.mpiexec()
    jlcmd = Base.julia_cmd()
    pdir = pkgdir(Peridynamics)
    vtk_path = joinpath(root, "vtk")

    # Test 1: Legacy mode - file I/O without result collection
    root_post_mpi = joinpath(root, "post_mpi")
    mkpath(root_post_mpi)

    mpi_cmd_legacy = """
    using Peridynamics
    files = "$(vtk_path)"
    process_each_export(files) do result0, result, file_id
        filename = string("max_displacement_", file_id, ".txt")
        open(joinpath("$(root_post_mpi)", filename), "w+") do io
            maxdisp = maximum(result[:displacement][1,:])
            msg = string("maximum displacement x: ", maxdisp)
            write(io, msg)
        end
        return nothing
    end
    """
    run(`$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd_legacy)`)
    @test isfile(joinpath(root_post_mpi, "max_displacement_1.txt"))
    @test contains(read(joinpath(root_post_mpi, "max_displacement_1.txt"), String),
                   "maximum displacement x: 0.0")
    @test isfile(joinpath(root_post_mpi, "max_displacement_6.txt"))

    # Test 2: MPI with barrier
    counter_file = joinpath(root, "counter.txt")
    rm(counter_file; force=true)

    mpi_cmd_barrier = """
    using Peridynamics
    files = "$(vtk_path)"
    counter_file = "$(counter_file)"

    # Test serial with barrier - all ranks should reach the barrier
    process_each_export(files; serial=true, barrier=true) do r0, r, id
        if mpi_isroot()
            open(counter_file, "a") do io
                println(io, "processed: \$id")
            end
        end
    end

    # If we get here, barrier worked correctly
    if mpi_isroot()
        open(counter_file, "a") do io
            println(io, "all ranks synchronized")
        end
    end
    """
    @test success(`$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd_barrier)`)
    @test isfile(counter_file)
    content = read(counter_file, String)
    @test contains(content, "processed: 1")
    @test contains(content, "all ranks synchronized")

    # Test 3: MPI parallel mode with result collection (serial=false)
    results_dir_parallel = joinpath(root, "mpi_results_parallel")
    mkpath(results_dir_parallel)

    mpi_cmd_parallel = """
    using Peridynamics, Test
    files = "$(vtk_path)"
    results_dir = "$(results_dir_parallel)"

    # Test result collection with MPI parallel mode (serial=false)
    default_value = (; max_disp=NaN, min_disp=NaN, file_id=0)
    results = process_each_export(files, default_value; serial=false) do r0, r, id
        max_disp = maximum(r[:displacement])
        min_disp = minimum(r[:displacement])
        return (; max_disp, min_disp, file_id=id)
    end

    # All ranks should have the same complete results
    @test length(results) == 6  # reference + 5 time steps
    @test all(results[i].file_id == i for i in 1:6)

    # Write results from each rank to verify they're identical
    rank = Peridynamics.mpi_rank()
    output_file = joinpath(results_dir, "parallel_rank_\$(rank).txt")
    open(output_file, "w") do io
        for r in results
            println(io, "file_id=\$(r.file_id), max=\$(r.max_disp), min=\$(r.min_disp)")
        end
    end
    """
    run_result = run(`$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd_parallel)`)
    @test run_result.exitcode == 0
    @test isfile(joinpath(results_dir_parallel, "parallel_rank_0.txt"))
    @test isfile(joinpath(results_dir_parallel, "parallel_rank_1.txt"))
    results_rank0 = read(joinpath(results_dir_parallel, "parallel_rank_0.txt"), String)
    results_rank1 = read(joinpath(results_dir_parallel, "parallel_rank_1.txt"), String)
    @test results_rank0 == results_rank1
    @test contains(results_rank0, "file_id=1")
    @test contains(results_rank0, "file_id=6")

    # Test 4: MPI serial mode with broadcast (serial=true)
    results_dir_serial = joinpath(root, "mpi_results_serial")
    mkpath(results_dir_serial)

    mpi_cmd_serial = """
    using Peridynamics, Test
    files = "$(vtk_path)"
    results_dir = "$(results_dir_serial)"

    # Test result collection with MPI serial mode (serial=true) - tests broadcast
    default_value = (; max_disp=NaN, avg_disp=NaN, file_id=0)
    results = process_each_export(files, default_value; serial=true) do r0, r, id
        max_disp = maximum(r[:displacement])
        avg_disp = sum(r[:displacement]) / length(r[:displacement])
        return (; max_disp, avg_disp, file_id=id)
    end

    # All ranks should have the same complete results (via broadcast)
    @test length(results) == 6
    @test all(results[i].file_id == i for i in 1:6)

    # Write results from each rank to verify broadcast worked
    rank = Peridynamics.mpi_rank()
    output_file = joinpath(results_dir, "serial_rank_\$(rank).txt")
    open(output_file, "w") do io
        for r in results
            println(io, "file_id=\$(r.file_id), max=\$(r.max_disp), avg=\$(r.avg_disp)")
        end
    end
    """
    run_result = run(`$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd_serial)`)
    @test run_result.exitcode == 0
    @test isfile(joinpath(results_dir_serial, "serial_rank_0.txt"))
    @test isfile(joinpath(results_dir_serial, "serial_rank_1.txt"))
    results_rank0 = read(joinpath(results_dir_serial, "serial_rank_0.txt"), String)
    results_rank1 = read(joinpath(results_dir_serial, "serial_rank_1.txt"), String)
    @test results_rank0 == results_rank1
    @test contains(results_rank0, "file_id=1")
    @test contains(results_rank0, "file_id=6")

    # Test 5: Non-bitstype error with MPI
    mpi_cmd_error = """
    using Peridynamics, Test
    files = "$(vtk_path)"

    # Test that non-bitstype default_value throws error with MPI
    default_value = (; name="test")  # String is not a bitstype
    @test_throws ArgumentError process_each_export(files, default_value) do r0, r, id
        return (; name="result")
    end
    """
    @test success(`$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd_error)`)
end

@testitem "process_each_export with result collection - error handling" begin
    root = mktempdir()
    l, Δx = 1.0, 1 / 4
    pos, vol = uniform_box(l, l, l, Δx)
    b1 = Body(BBMaterial(), pos, vol)
    material!(b1; horizon=3.015Δx, E=2.1e5, rho=8e-6)
    velocity_ic!(b1, :all_points, :x, 1.0)
    vv = VelocityVerlet(steps=2)
    job = Job(b1, vv; path=root, freq=1)
    submit(job)

    # Test that errors during processing use default_value
    default_value = (; value=999.0)
    results = process_each_export(job, default_value) do r0, r, id
        if id == 2
            error("Intentional error for testing")
        end
        return (; value=Float64(id))
    end

    @test results[1].value == 1.0
    @test results[2].value == 999.0  # Should use default_value due to error
    @test results[3].value == 3.0

    # Verify error log was created (in vtk/process_errors)
    error_dir = joinpath(root, "vtk", "vtk_process_errors")
    @test isdir(error_dir)

    error_logs = filter(f -> contains(f, "proc_error_file"), readdir(error_dir))
    @test length(error_logs) == 1

    # Verify log content
    log_content = read(joinpath(error_dir, error_logs[1]), String)
    @test contains(log_content, "Intentional error for testing")
    @test contains(log_content, "File ID: 2")
end

@testitem "process_each_export error logging with serial mode" begin
    root = mktempdir()
    l, Δx = 1.0, 1 / 4
    pos, vol = uniform_box(l, l, l, Δx)
    b1 = Body(BBMaterial(), pos, vol)
    material!(b1; horizon=3.015Δx, E=2.1e5, rho=8e-6)
    velocity_ic!(b1, :all_points, :x, 1.0)
    vv = VelocityVerlet(steps=5)
    job = Job(b1, vv; path=root, freq=1)
    submit(job)

    # Test multiple errors in serial mode
    default_value = (; max_disp=0.0)
    results = process_each_export(job, default_value; serial=true) do r0, r, id
        if id == 2 || id == 4
            error("Error at timestep $id")
        end
        return (; max_disp=Float64(id) * 0.1)
    end

    @test results[1].max_disp ≈ 0.1
    @test results[2].max_disp == 0.0  # Error -> default
    @test results[3].max_disp ≈ 0.3
    @test results[4].max_disp == 0.0  # Error -> default
    @test results[5].max_disp ≈ 0.5
    @test results[6].max_disp ≈ 0.6

    # Verify error logs were created for both errors (in vtk/process_errors)
    error_dir = joinpath(root, "vtk", "vtk_process_errors")
    @test isdir(error_dir)

    error_logs = filter(f -> contains(f, "proc_error_file"), readdir(error_dir))
    @test length(error_logs) == 2

    # Verify both error logs contain correct information
    log_contents = [read(joinpath(error_dir, log), String) for log in error_logs]
    @test any(contains(content, "File ID: 2") for content in log_contents)
    @test any(contains(content, "File ID: 4") for content in log_contents)
end

@testitem "process_each_export error logging without result collection" begin
    root = mktempdir()
    l, Δx = 1.0, 1 / 4
    pos, vol = uniform_box(l, l, l, Δx)
    b1 = Body(BBMaterial(), pos, vol)
    material!(b1; horizon=3.015Δx, E=2.1e5, rho=8e-6)
    velocity_ic!(b1, :all_points, :x, 1.0)
    vv = VelocityVerlet(steps=3)
    job = Job(b1, vv; path=root, freq=1)
    submit(job)

    # Test error logging in legacy mode (no result collection)
    counter = Ref(0)
    result = process_each_export(job; serial=true) do r0, r, id
        counter[] += 1
        if id == 3
            error("Legacy mode error test")
        end
        return nothing
    end

    @test result === nothing
    @test counter[] == 4  # All files processed despite error

    # Verify error log was created (in vtk/process_errors)
    error_dir = joinpath(root, "vtk", "vtk_process_errors")
    @test isdir(error_dir)

    error_logs = filter(f -> contains(f, "proc_error_file"), readdir(error_dir))
    @test length(error_logs) == 1

    log_content = read(joinpath(error_dir, error_logs[1]), String)
    @test contains(log_content, "Legacy mode error test")
    @test contains(log_content, "File ID: 3")
end



@testitem "process_each_export backwards compatibility" begin
    root = mktempdir()
    l, Δx = 1.0, 1 / 4
    pos, vol = uniform_box(l, l, l, Δx)
    b1 = Body(BBMaterial(), pos, vol)
    material!(b1; horizon=3.015Δx, E=2.1e5, rho=8e-6)
    velocity_ic!(b1, :all_points, :x, 1.0)
    vv = VelocityVerlet(steps=1)
    job = Job(b1, vv; path=root, freq=1)
    submit(job)

    # Test that old behavior still works (returns nothing)
    counter = Ref(0)
    result = process_each_export(job; serial=true) do r0, r, id
        counter[] += 1
        return nothing
    end

    @test result === nothing
    @test counter[] == 2  # reference + 1 time step

    # Test with explicit nothing default_value
    counter[] = 0
    result2 = process_each_export(job, nothing; serial=true) do r0, r, id
        counter[] += 1
        return nothing
    end

    @test result2 === nothing
    @test counter[] == 2
end

@testitem "process_each_export with different result types" begin
    root = mktempdir()
    l, Δx = 1.0, 1 / 4
    pos, vol = uniform_box(l, l, l, Δx)
    b1 = Body(BBMaterial(), pos, vol)
    material!(b1; horizon=3.015Δx, E=2.1e5, rho=8e-6)
    velocity_ic!(b1, :all_points, :x, 1.0)
    vv = VelocityVerlet(steps=2)
    job = Job(b1, vv; path=root, freq=1)
    submit(job)

    # Test with Int result
    results_int = process_each_export(job, 0) do r0, r, id
        return id
    end
    @test results_int isa Vector{Int}
    @test results_int == [1, 2, 3]

    # Test with Float64 result
    results_float = process_each_export(job, 0.0) do r0, r, id
        return Float64(id) * 1.5
    end
    @test results_float isa Vector{Float64}
    @test results_float ≈ [1.5, 3.0, 4.5]

    # Test with Bool result
    results_bool = process_each_export(job, false) do r0, r, id
        return id > 2
    end
    @test results_bool isa Vector{Bool}
    @test results_bool == [false, false, true]

    # Test with complex NamedTuple
    default = (; a=0.0, b=0, c=false)
    results_complex = process_each_export(job, default) do r0, r, id
        return (; a=Float64(id), b=id*2, c=(id == 2))
    end
    @test results_complex isa Vector{NamedTuple{(:a, :b, :c), Tuple{Float64, Int, Bool}}}
    @test results_complex[1] == (; a=1.0, b=2, c=false)
    @test results_complex[2] == (; a=2.0, b=4, c=true)
    @test results_complex[3] == (; a=3.0, b=6, c=false)
end
