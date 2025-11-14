@testitem "process_each_export - serial and threads" begin
    # Single simulation for all serial and threading tests
    Peridynamics.MPI_RUN[] = false

    root = mktempdir()
    l, Δx = 1.0, 1 / 4
    pos, vol = uniform_box(l, l, l, Δx)
    body = Body(BBMaterial(), pos, vol)
    material!(body; horizon=3.015Δx, E=2.1e5, rho=8e-6)
    point_set!(y -> y > l / 2 - Δx, body, :set_top)
    point_set!(y -> y < -l / 2 + Δx, body, :set_bottom)
    velocity_bc!(t -> 30, body, :set_top, :y)
    velocity_bc!(t -> -30, body, :set_bottom, :y)
    vv = VelocityVerlet(steps=5)
    job = Job(body, vv; path=root, freq=1)

    # Test argument validation before submission
    @test_throws ArgumentError process_each_export((r0, r, id) -> nothing, job)

    submit(job)

    # Test 1: Legacy mode with file I/O (threads backend)
    root_post_threads = joinpath(root, "post_threads")
    mkpath(root_post_threads)
    process_each_export(job) do r0, r, id
        filename = "max_disp_$(id).txt"
        open(joinpath(root_post_threads, filename), "w") do io
            write(io, "max_disp: $(maximum(r[:displacement][1,:]))")
        end
        return nothing
    end
    @test isfile(joinpath(root_post_threads, "max_disp_1.txt"))
    @test contains(read(joinpath(root_post_threads, "max_disp_1.txt"), String), "max_disp: 0.0")
    @test isfile(joinpath(root_post_threads, "max_disp_6.txt"))

    # Test 2: Legacy mode (serial backend)
    root_post_serial = joinpath(root, "post_serial")
    mkpath(root_post_serial)
    process_each_export(job; serial=true) do r0, r, id
        filename = "max_disp_$(id).txt"
        open(joinpath(root_post_serial, filename), "w") do io
            write(io, "max_disp: $(maximum(r[:displacement][1,:]))")
        end
        return nothing
    end
    @test isfile(joinpath(root_post_serial, "max_disp_1.txt"))
    @test contains(read(joinpath(root_post_serial, "max_disp_1.txt"), String), "max_disp: 0.0")
    @test isfile(joinpath(root_post_serial, "max_disp_6.txt"))

    # Test 3: Barrier parameter (non-MPI context)
    counter = Ref(0)
    process_each_export(job; serial=true, barrier=false) do r0, r, id
        counter[] += 1
    end
    @test counter[] == 6  # reference + 5 time steps

    counter[] = 0
    process_each_export(job; serial=true, barrier=true) do r0, r, id
        counter[] += 1
    end
    @test counter[] == 6

    # Test 4: Result collection with NamedTuples (threads backend)
    default_nt = (; max_disp=NaN, file_id=0)
    results_threads = process_each_export(job, default_nt; serial=false) do r0, r, id
        return (; max_disp=maximum(r[:displacement]), file_id=id)
    end
    @test results_threads isa Vector{NamedTuple{(:max_disp, :file_id), Tuple{Float64, Int}}}
    @test length(results_threads) == 6
    @test all(results_threads[i].file_id == i for i in 1:6)
    @test results_threads[1].max_disp ≈ 0.0
    @test results_threads[end].max_disp > 0.0

    # Test 5: Result collection (serial backend)
    default_avg = (; avg_ux=NaN, avg_uy=NaN)
    results_serial = process_each_export(job, default_avg; serial=true) do r0, r, id
        ux = @view r[:displacement][1, :]
        uy = @view r[:displacement][2, :]
        return (; avg_ux=sum(ux)/length(ux), avg_uy=sum(uy)/length(uy))
    end
    @test results_serial isa Vector{NamedTuple{(:avg_ux, :avg_uy), Tuple{Float64, Float64}}}
    @test length(results_serial) == 6
    @test all(.!isnan.(getfield.(results_serial, :avg_ux)))
    @test results_serial[1].avg_ux ≈ 0.0
    @test results_serial[1].avg_uy ≈ 0.0

    # Test 6: Different result types
    results_int = process_each_export(job, 0) do r0, r, id
        return id
    end
    @test results_int isa Vector{Int}
    @test results_int == [1, 2, 3, 4, 5, 6]

    results_float = process_each_export(job, 0.0) do r0, r, id
        return Float64(id) * 2.5
    end
    @test results_float isa Vector{Float64}
    @test results_float == [2.5, 5.0, 7.5, 10.0, 12.5, 15.0]

    results_bool = process_each_export(job, false) do r0, r, id
        return id > 3
    end
    @test results_bool isa Vector{Bool}
    @test results_bool == [false, false, false, true, true, true]

    default_complex = (; a=0.0, b=0, c=false)
    results_complex = process_each_export(job, default_complex) do r0, r, id
        return (; a=Float64(id), b=id*2, c=(id == 3))
    end
    @test results_complex isa Vector{NamedTuple{(:a, :b, :c), Tuple{Float64, Int, Bool}}}
    @test results_complex[3] == (; a=3.0, b=6, c=true)

    # Test 7: Error handling with result collection
    default_err = (; value=999.0)
    results_err = process_each_export(job, default_err) do r0, r, id
        if id == 2
            error("Intentional error for testing")
        end
        return (; value=Float64(id))
    end
    @test results_err[1].value == 1.0
    @test results_err[2].value == 999.0  # Uses default due to error
    @test results_err[3].value == 3.0

    error_dir = joinpath(root, "vtk", "vtk_process_errors")
    @test isdir(error_dir)
    error_logs = filter(f -> contains(f, "proc_error_file"), readdir(error_dir))
    @test length(error_logs) == 1
    log_content = read(joinpath(error_dir, error_logs[1]), String)
    @test contains(log_content, "Intentional error for testing")
    @test contains(log_content, "File ID: 2")

    # Test 8: Multiple errors (serial mode)
    default_multi = (; val=0.0)
    results_multi = process_each_export(job, default_multi; serial=true) do r0, r, id
        if id == 3 || id == 5
            error("Error at timestep $id")
        end
        return (; val=Float64(id) * 0.1)
    end
    @test results_multi[2].val ≈ 0.2
    @test results_multi[3].val == 0.0  # Error -> default
    @test results_multi[4].val ≈ 0.4
    @test results_multi[5].val == 0.0  # Error -> default

    error_logs = filter(f -> contains(f, "proc_error_file"), readdir(error_dir))
    @test length(error_logs) == 3  # 1 from previous test + 2 new

    # Test 9: Error logging without result collection
    counter[] = 0
    process_each_export(job; serial=true) do r0, r, id
        counter[] += 1
        if id == 4
            error("Legacy mode error test")
        end
        return nothing
    end
    @test counter[] == 6  # All files processed despite error

    error_logs = filter(f -> contains(f, "proc_error_file"), readdir(error_dir))
    @test length(error_logs) == 4  # 3 from previous + 1 new

    # Test 10: Backwards compatibility
    result_nothing = process_each_export(job) do r0, r, id
        return nothing
    end
    @test result_nothing === nothing

    result_nothing2 = process_each_export(job, nothing) do r0, r, id
        return nothing
    end
    @test result_nothing2 === nothing
end

@testitem "find_vtk_files" begin
    root = mktempdir()
    vtk_files_unsorted = ["timestep_123.pvtu", "abcd_timestep_000005.pvtu",
                          "timestep_02.pvtu", "_1.pvtu"]
    for file in vtk_files_unsorted
        open(joinpath(root, file), "w") do io
            write(io, "")
        end
    end

    vtk_files = Peridynamics.find_vtk_files(root)
    @test basename.(vtk_files) == ["_1.pvtu", "timestep_02.pvtu",
                                   "abcd_timestep_000005.pvtu", "timestep_123.pvtu"]
end

@testitem "process_each_export - MPI tests" tags=[:mpi] begin
    # Single simulation for all MPI tests
    root = mktempdir()
    l, Δx = 1.0, 1 / 4
    pos, vol = uniform_box(l, l, l, Δx)
    body = Body(BBMaterial(), pos, vol)
    material!(body; horizon=3.015Δx, E=2.1e5, rho=8e-6)
    point_set!(y -> y > l / 2 - Δx, body, :set_top)
    point_set!(y -> y < -l / 2 + Δx, body, :set_bottom)
    velocity_bc!(t -> 30, body, :set_top, :y)
    velocity_bc!(t -> -30, body, :set_bottom, :y)
    vv = VelocityVerlet(steps=5)
    job = Job(body, vv; path=root, freq=1)
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
    process_each_export(files) do r0, r, id
        filename = "max_disp_\$(id).txt"
        open(joinpath("$(root_post_mpi)", filename), "w") do io
            write(io, "max_disp: \$(maximum(r[:displacement][1,:]))")
        end
        return nothing
    end
    """
    run(`$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd_legacy)`)
    @test isfile(joinpath(root_post_mpi, "max_disp_1.txt"))
    @test contains(read(joinpath(root_post_mpi, "max_disp_1.txt"), String), "max_disp: 0.0")
    @test isfile(joinpath(root_post_mpi, "max_disp_6.txt"))

    # Test 2: MPI with barrier
    counter_file = joinpath(root, "counter.txt")
    mpi_cmd_barrier = """
    using Peridynamics
    files = "$(vtk_path)"
    counter_file = "$(counter_file)"
    process_each_export(files; serial=true, barrier=true) do r0, r, id
        if mpi_isroot()
            open(counter_file, "a") do io
                println(io, "processed: \$id")
            end
        end
    end
    if mpi_isroot()
        open(counter_file, "a") do io
            println(io, "all ranks synchronized")
        end
    end
    """
    @test success(`$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd_barrier)`)
    @test isfile(counter_file)
    @test contains(read(counter_file, String), "processed: 1")
    @test contains(read(counter_file, String), "all ranks synchronized")

    # Test 3: MPI parallel mode with result collection (serial=false)
    results_dir_parallel = joinpath(root, "mpi_results_parallel")
    mkpath(results_dir_parallel)
    mpi_cmd_parallel = """
    using Peridynamics, Test
    files = "$(vtk_path)"
    results_dir = "$(results_dir_parallel)"
    default_value = (; max_disp=NaN, min_disp=NaN, file_id=0)
    results = process_each_export(files, default_value; serial=false) do r0, r, id
        return (; max_disp=maximum(r[:displacement]), min_disp=minimum(r[:displacement]), file_id=id)
    end
    @test length(results) == 6
    @test all(results[i].file_id == i for i in 1:6)
    rank = Peridynamics.mpi_rank()
    open(joinpath(results_dir, "parallel_rank_\$(rank).txt"), "w") do io
        for r in results
            println(io, "file_id=\$(r.file_id), max=\$(r.max_disp), min=\$(r.min_disp)")
        end
    end
    """
    @test success(`$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd_parallel)`)
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
    default_value = (; max_disp=NaN, avg_disp=NaN, file_id=0)
    results = process_each_export(files, default_value; serial=true) do r0, r, id
        disp = r[:displacement]
        return (; max_disp=maximum(disp), avg_disp=sum(disp)/length(disp), file_id=id)
    end
    @test length(results) == 6
    @test all(results[i].file_id == i for i in 1:6)
    rank = Peridynamics.mpi_rank()
    open(joinpath(results_dir, "serial_rank_\$(rank).txt"), "w") do io
        for r in results
            println(io, "file_id=\$(r.file_id), max=\$(r.max_disp), avg=\$(r.avg_disp)")
        end
    end
    """
    @test success(`$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd_serial)`)
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
    default_value = (; name="test")  # String is not a bitstype
    @test_throws ArgumentError process_each_export(files, default_value) do r0, r, id
        return (; name="result")
    end
    """
    @test success(`$(mpiexec) -n 2 $(jlcmd) --project=$(pdir) -e $(mpi_cmd_error)`)
end
