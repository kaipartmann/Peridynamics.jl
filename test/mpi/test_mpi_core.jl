@testitem "MPI: core paths" tags=[:mpi] setup=[Fixtures] begin
    # One `mpiexec` run for every MPI code path that needs two ranks but no big simulation:
    # halo exchanges, `@mpiroot`, progress bars, the synchronized `NaNError`, `process_each_export`,
    # `Study` with `only_root` and the MPI timers. See `scripts/core.jl` for the checks; what
    # the script writes to disk is verified here.
    root = mktempdir()
    mkpath(joinpath(root, "post"))

    # a finished threads simulation whose VTK files the script post-processes
    vtk_root = joinpath(root, "sim")
    l, Δx = 1.0, 1 / 4
    pos, vol = uniform_box(l, l, l, Δx)
    body = Body(BBMaterial(), pos, vol)
    material!(body; horizon=3.015Δx, E=2.1e5, rho=8e-6)
    point_set!(y -> y > l / 2 - Δx, body, :set_top)
    point_set!(y -> y < -l / 2 + Δx, body, :set_bottom)
    velocity_bc!(t -> 30, body, :set_top, :y)
    velocity_bc!(t -> -30, body, :set_bottom, :y)
    submit(Job(body, VelocityVerlet(steps=5); path=vtk_root, freq=1); quiet=true)
    vtk_path = joinpath(vtk_root, "vtk")
    @test length(Peridynamics.find_vtk_files(vtk_path)) == 6

    @test Fixtures.run_mpi_script("core.jl", root, vtk_path; nranks=2)

    # NaNError: the message is in the logfile of the aborted job
    logfile = read(joinpath(root, "nan", "logfile.log"), String)
    Δt = 1e-7
    @test occursin("NaN values detected in force density field!\n  time:    $(2Δt)\n  step:    2\n",
                   logfile)

    # process_each_export: the files written by the ranks
    post = joinpath(root, "post")
    @test contains(read(joinpath(post, "max_disp_1.txt"), String), "max_disp: 0.0")
    @test isfile(joinpath(post, "max_disp_6.txt"))
    counter = read(joinpath(post, "counter.txt"), String)
    @test contains(counter, "processed: 1") && contains(counter, "processed: 6")
    @test contains(counter, "all ranks synchronized")
    @test !contains(counter, "processed: 7")
    for mode in ("parallel", "serial")
        rank0 = read(joinpath(post, "$(mode)_rank_0.txt"), String)
        rank1 = read(joinpath(post, "$(mode)_rank_1.txt"), String)
        @test rank0 == rank1
        @test contains(rank0, "file_id=1") && contains(rank0, "file_id=6")
    end

    # Study with only_root: only the root rank processed
    marker = read(joinpath(root, "processing_marker.txt"), String)
    @test contains(marker, "Job 1 processed on rank 0")
    @test contains(marker, "Job 2 processed on rank 0")
    @test !contains(marker, "rank 1")

    # timers: one log per rank
    @test isfile(joinpath(root, "timers", "timers_rank_0.log"))
    @test isfile(joinpath(root, "timers", "timers_rank_1.log"))
end
