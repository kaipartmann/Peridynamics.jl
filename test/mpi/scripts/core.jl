# The MPI code paths that need more than one rank, checked in a single `mpiexec` run:
#     mpiexec -n 2 julia --project=<Peridynamics> test/mpi/scripts/core.jl <root> <vtk_path>
# `<root>` is a writable scratch directory, `<vtk_path>` a directory with the VTK output of a
# finished simulation (for `process_each_export`). Every check is a `@testset`, so a failure
# names the check and the script exits non-zero at the end of that testset. What the parent test
# item verifies from the outside (log files, marker files) is written below `<root>`.
#
# One process instead of one per check, because with coverage on every fresh Julia process pays a
# full instrumented compilation of the package; that is the whole cost of these tests.
using Peridynamics, Test

const ROOT = ARGS[1]
const VTK_PATH = ARGS[2]
const RANK = Peridynamics.mpi_rank()

@assert Peridynamics.mpi_run()
@assert Peridynamics.mpi_nranks() == 2

# the 4-point body of the halo exchange checks: every point sees every other point
function tetra4()
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(CMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1, nu=0.25, Gc=1)
    return body
end

const TESTPOS = [0.1 0.2 0.3 0.4
                 0.1 0.2 0.3 0.4
                 0.1 0.2 0.3 0.4]
const TESTBINT = [0.11 0.21 0.31 0.41
                  0.12 0.22 0.32 0.42
                  0.13 0.23 0.33 0.43]

@testset "exchange_loc_to_halo!" begin
    dh = Peridynamics.mpi_data_handler(tetra4(), VelocityVerlet(steps=10))
    if RANK == 1
        dh.chunk.storage.position .= TESTPOS
        dh.chunk.storage.b_int .= TESTBINT
    end
    Peridynamics.exchange_loc_to_halo!(dh)
    Peridynamics.exchange_loc_to_halo!(chunk -> chunk.storage.b_int, dh)
    if RANK == 0
        @test dh.chunk.storage.position[:, 3:4] ≈ TESTPOS[:, 1:2]
        @test dh.chunk.storage.b_int[:, 3:4] ≈ TESTBINT[:, 1:2]
    end
end

@testset "exchange_halo_to_loc!" begin
    body = tetra4()
    dh = Peridynamics.mpi_data_handler(body, VelocityVerlet(steps=10))
    if RANK == 1
        dh.chunk.storage.position .= TESTPOS
        dh.chunk.storage.b_int .= TESTBINT
    end
    if RANK == 0 # the destination chunk is untouched before the exchange
        @test iszero(dh.chunk.storage.b_int)
        @test dh.chunk.storage.position ≈ body.position
    end
    Peridynamics.exchange_halo_to_loc!(dh)
    Peridynamics.exchange_halo_to_loc!(chunk -> chunk.storage.position, dh)
    if RANK == 0
        @test dh.chunk.storage.b_int[:, 1:2] ≈ TESTBINT[:, 3:4]
        @test dh.chunk.storage.position[:, 1:2] ≈ body.position[:, 1:2] + TESTPOS[:, 3:4]
    end
end

@testset "@mpiroot and progress bars" begin
    counter = Ref(0)
    @mpiroot counter[] += 1
    @test counter[] == (RANK == 0 ? 1 : 0)
    @mpiroot :wait counter[] += 1
    @test counter[] == (RANK == 0 ? 2 : 0)
    @test mpi_isroot() == (RANK == 0)
    # progress bars are off with MPI unless enabled by hand, which warns on every init
    @test !Peridynamics.mpi_progress_bars()
    enable_mpi_progress_bars!()
    @test Peridynamics.mpi_progress_bars()
    @test_logs (:warn, r"progress bar settings overwritten") Peridynamics.set_progress_bars!()
    reset_mpi_progress_bars!()
    @test !Peridynamics.mpi_progress_bars()
    Peridynamics.set_progress_bars!()
    @test !Peridynamics.progress_bars()
end

@testset "submit: NaNError is rethrown on all ranks" begin
    # the velocity jumps to 1e40 in step 2, so the force density contains NaNs in step 2; the
    # NaNError is synchronized across the ranks and rethrown on every rank
    path = joinpath(ROOT, "nan")
    Δt = 1e-7
    pos = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    body = Body(BBMaterial(), pos, [1.0, 1.0])
    material!(body, horizon=1.5, rho=8000, E=210e9)
    velocity_bc!(t -> t < 1.9Δt ? 1.0 : 1e40, body, :all_points, :x)
    job = Job(body, VelocityVerlet(steps=5, stepsize=Δt); path)
    @test_throws Peridynamics.NaNError(2Δt, 2) submit(job; quiet=true)
end

@testset "process_each_export" begin
    files = VTK_PATH
    post = joinpath(ROOT, "post")
    # file output without result collection
    process_each_export(files) do r0, r, id
        open(joinpath(post, "max_disp_$(id).txt"), "w") do io
            write(io, "max_disp: $(maximum(r[:displacement][1, :]))")
        end
        return nothing
    end
    # serial with barrier: only root processes, and all ranks wait for it
    counter_file = joinpath(post, "counter.txt")
    process_each_export(files; serial=true, barrier=true) do r0, r, id
        if mpi_isroot()
            open(counter_file, "a") do io
                println(io, "processed: $id")
            end
        end
    end
    @mpiroot open(counter_file, "a") do io
        println(io, "all ranks synchronized")
    end
    # parallel mode with result collection: every rank ends up with all results
    default_value = (; max_disp=NaN, min_disp=NaN, file_id=0)
    results = process_each_export(files, default_value; serial=false) do r0, r, id
        return (; max_disp=maximum(r[:displacement]), min_disp=minimum(r[:displacement]),
                  file_id=id)
    end
    @test length(results) == 6
    @test all(results[i].file_id == i for i in 1:6)
    open(joinpath(post, "parallel_rank_$(RANK).txt"), "w") do io
        for r in results
            println(io, "file_id=$(r.file_id), max=$(r.max_disp), min=$(r.min_disp)")
        end
    end
    # serial mode with broadcast of the results
    default_value2 = (; max_disp=NaN, avg_disp=NaN, file_id=0)
    results2 = process_each_export(files, default_value2; serial=true) do r0, r, id
        disp = r[:displacement]
        return (; max_disp=maximum(disp), avg_disp=sum(disp) / length(disp), file_id=id)
    end
    @test length(results2) == 6
    @test all(results2[i].file_id == i for i in 1:6)
    open(joinpath(post, "serial_rank_$(RANK).txt"), "w") do io
        for r in results2
            println(io, "file_id=$(r.file_id), max=$(r.max_disp), avg=$(r.avg_disp)")
        end
    end
    # results must be bits types to be communicated
    @test_throws ArgumentError process_each_export(files, (; name="test")) do r0, r, id
        return (; name="result")
    end
end

@testset "Study: submit! and process_each_job with only_root" begin
    function create_job(setup::NamedTuple, root::String)
        position = zeros(3, 10)
        position[1, :] .= 0.0:9.0
        body = Body(BBMaterial(), position, ones(10))
        material!(body, horizon=1.5, E=1.0, rho=1, Gc=1)
        velocity_ic!(body, :all_points, :x, 1.0)
        return Job(body, VelocityVerlet(steps=5); path=joinpath(root, "sim_$(setup.id)"), freq=5)
    end
    study_root = joinpath(ROOT, "study")
    study = Study(create_job, [(; id=1), (; id=2)]; root=study_root)
    submit!(study; quiet=true)
    @test all(study.sim_success)
    marker_file = joinpath(ROOT, "processing_marker.txt")
    results = process_each_job(study, (; id=0, rank=-1); only_root=true) do job, setup
        if mpi_isroot() # this should only run on root
            open(marker_file, "a") do io
                println(io, "Job $(setup.id) processed on rank $(RANK)")
            end
        end
        return (; id=setup.id, rank=RANK)
    end
    if mpi_isroot()
        @test isfile(marker_file)
        content = read(marker_file, String)
        @test contains(content, "Job 1 processed on rank 0")
        @test contains(content, "Job 2 processed on rank 0")
        @test !contains(content, "rank 1")
    end
end

# The timers are enabled by evaluating a method definition into the package, so the call that
# uses them must be a separate top-level statement (a later world age) — not inside the same
# `@testset` block as `enable_mpi_timers!()`.
disable_mpi_timers!()
enable_mpi_timers!()
@testset "MPI timers" begin
    path = joinpath(ROOT, "timers")
    body = tetra4()
    velocity_ic!(body, :all_points, :x, 1.0)
    job = Job(body, VelocityVerlet(steps=10); path, freq=10)
    submit(job; quiet=true)
    @mpiroot :wait nothing
    @test isfile(joinpath(path, "timers_rank_0.log"))
    @test isfile(joinpath(path, "timers_rank_1.log"))
end
disable_mpi_timers!()
