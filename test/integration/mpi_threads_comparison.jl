@testitem "MPI-Threads comparison BBMaterial{NoCorrection}" begin
    root = joinpath(@__DIR__, "temp_mpi_threads_comparison_bbnc")
    path_threads = joinpath(root, "results_threads")
    path_threads_vtk = joinpath(path_threads, "vtk")
    path_mpi = joinpath(root, "results_mpi")
    path_mpi_vtk = joinpath(path_mpi, "vtk")

    function sim_bb(N::Int, path::String)
        l, Δx, δ, a = 1.0, 1/N, 3.015/N, 0.5
        pos, vol = uniform_box(l, l, 0.1l, Δx)
        ids = sortperm(pos[2,:])
        b = Body(BBMaterial(), pos[:, ids], vol[ids])
        material!(b; horizon=3.015Δx, E=2.1e5, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
        point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
        precrack!(b, :set_a, :set_b)
        point_set!(p -> p[2] > l/2-Δx, b, :set_top)
        point_set!(p -> p[2] < -l/2+Δx, b, :set_bottom)
        velocity_bc!(t -> -30, b, :set_bottom, :y)
        velocity_bc!(t -> 30, b, :set_top, :y)
        vv = VelocityVerlet(steps=100)
        job = Job(b, vv; path=path, freq=50)
        submit(job)
        return nothing
    end
    sim_bb(30, path_threads)

    mpi_cmd = """
    using Peridynamics
    function sim_bb(N::Int, path::String)
        l, Δx, δ, a = 1.0, 1/N, 3.015/N, 0.5
        pos, vol = uniform_box(l, l, 0.1l, Δx)
        ids = sortperm(pos[2,:])
        b = Body(BBMaterial(), pos[:, ids], vol[ids])
        material!(b; horizon=3.015Δx, E=2.1e5, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
        point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
        precrack!(b, :set_a, :set_b)
        point_set!(p -> p[2] > l/2-Δx, b, :set_top)
        point_set!(p -> p[2] < -l/2+Δx, b, :set_bottom)
        velocity_bc!(t -> -30, b, :set_bottom, :y)
        velocity_bc!(t -> 30, b, :set_top, :y)
        vv = VelocityVerlet(steps=100)
        job = Job(b, vv; path=path, freq=50)
        submit(job)
        return nothing
    end
    sim_bb(30, "$path_mpi")
    """
    run(`$(Peridynamics.MPI.mpiexec()) -n 2 $(Base.julia_cmd()) --project -e $(mpi_cmd)`)

    @test isdir(path_threads_vtk)
    @test isdir(path_mpi_vtk)
    vtk_files_threads = Peridynamics.find_vtk_files(path_threads_vtk)
    vtk_files_mpi = Peridynamics.find_vtk_files(path_mpi_vtk)
    @test length(vtk_files_mpi) == length(vtk_files_threads) == 3
    for i in eachindex(vtk_files_threads, vtk_files_mpi)
        res_threads = read_vtk(vtk_files_threads[i])
        res_mpi = read_vtk(vtk_files_mpi[i])
        for key in keys(res_threads)
            @test res_threads[key] ≈ res_mpi[key]
        end
    end

    rm(root; recursive=true, force=true)
end

@testitem "MPI-Threads comparison BBMaterial{EnergySurfaceCorrection}" begin
    root = joinpath(@__DIR__, "temp_mpi_threads_comparison_bbesc")
    path_threads = joinpath(root, "results_threads")
    path_threads_vtk = joinpath(path_threads, "vtk")
    path_mpi = joinpath(root, "results_mpi")
    path_mpi_vtk = joinpath(path_mpi, "vtk")

    function sim_bb(N::Int, path::String)
        l, Δx, δ, a = 1.0, 1/N, 3.015/N, 0.5
        pos, vol = uniform_box(l, l, 0.1l, Δx)
        ids = sortperm(pos[2,:])
        b = Body(BBMaterial{EnergySurfaceCorrection}(), pos[:, ids], vol[ids])
        material!(b; horizon=3.015Δx, E=2.1e5, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
        point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
        precrack!(b, :set_a, :set_b)
        point_set!(p -> p[2] > l/2-Δx, b, :set_top)
        point_set!(p -> p[2] < -l/2+Δx, b, :set_bottom)
        velocity_bc!(t -> -30, b, :set_bottom, :y)
        velocity_bc!(t -> 30, b, :set_top, :y)
        vv = VelocityVerlet(steps=100)
        job = Job(b, vv; path=path, freq=50)
        submit(job)
        return nothing
    end
    sim_bb(30, path_threads)

    mpi_cmd = """
    using Peridynamics
    function sim_bb(N::Int, path::String)
        l, Δx, δ, a = 1.0, 1/N, 3.015/N, 0.5
        pos, vol = uniform_box(l, l, 0.1l, Δx)
        ids = sortperm(pos[2,:])
        b = Body(BBMaterial{EnergySurfaceCorrection}(), pos[:, ids], vol[ids])
        material!(b; horizon=3.015Δx, E=2.1e5, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
        point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
        precrack!(b, :set_a, :set_b)
        point_set!(p -> p[2] > l/2-Δx, b, :set_top)
        point_set!(p -> p[2] < -l/2+Δx, b, :set_bottom)
        velocity_bc!(t -> -30, b, :set_bottom, :y)
        velocity_bc!(t -> 30, b, :set_top, :y)
        vv = VelocityVerlet(steps=100)
        job = Job(b, vv; path=path, freq=50)
        submit(job)
        return nothing
    end
    sim_bb(30, "$path_mpi")
    """
    run(`$(Peridynamics.MPI.mpiexec()) -n 2 $(Base.julia_cmd()) --project -e $(mpi_cmd)`)

    @test isdir(path_threads_vtk)
    @test isdir(path_mpi_vtk)
    vtk_files_threads = Peridynamics.find_vtk_files(path_threads_vtk)
    vtk_files_mpi = Peridynamics.find_vtk_files(path_mpi_vtk)
    @test length(vtk_files_mpi) == length(vtk_files_threads) == 3
    for i in eachindex(vtk_files_threads, vtk_files_mpi)
        res_threads = read_vtk(vtk_files_threads[i])
        res_mpi = read_vtk(vtk_files_mpi[i])
        for key in keys(res_threads)
            @test res_threads[key] ≈ res_mpi[key]
        end
    end

    rm(root; recursive=true, force=true)
end

@testitem "MPI-Threads comparison OSBMaterial" begin
    root = joinpath(@__DIR__, "temp_mpi_threads_comparison_osb")
    path_threads = joinpath(root, "results_threads")
    path_threads_vtk = joinpath(path_threads, "vtk")
    path_mpi = joinpath(root, "results_mpi")
    path_mpi_vtk = joinpath(path_mpi, "vtk")

    function sim_osb(N::Int, path::String)
        l, Δx, δ, a = 1.0, 1/N, 3.015/N, 0.5
        pos, vol = uniform_box(l, l, 0.1l, Δx)
        ids = sortperm(pos[2,:])
        b = Body(OSBMaterial(), pos[:, ids], vol[ids])
        material!(b; horizon=3.015Δx, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
        point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
        precrack!(b, :set_a, :set_b)
        point_set!(p -> p[2] > l/2-Δx, b, :set_top)
        point_set!(p -> p[2] < -l/2+Δx, b, :set_bottom)
        velocity_bc!(t -> -30, b, :set_bottom, :y)
        velocity_bc!(t -> 30, b, :set_top, :y)
        vv = VelocityVerlet(steps=100)
        job = Job(b, vv; path=path, freq=50)
        submit(job)
        return nothing
    end
    sim_osb(30, path_threads)

    mpi_cmd = """
    using Peridynamics
    function sim_osb(N::Int, path::String)
        l, Δx, δ, a = 1.0, 1/N, 3.015/N, 0.5
        pos, vol = uniform_box(l, l, 0.1l, Δx)
        ids = sortperm(pos[2,:])
        b = Body(OSBMaterial(), pos[:, ids], vol[ids])
        material!(b; horizon=3.015Δx, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
        point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
        precrack!(b, :set_a, :set_b)
        point_set!(p -> p[2] > l/2-Δx, b, :set_top)
        point_set!(p -> p[2] < -l/2+Δx, b, :set_bottom)
        velocity_bc!(t -> -30, b, :set_bottom, :y)
        velocity_bc!(t -> 30, b, :set_top, :y)
        vv = VelocityVerlet(steps=100)
        job = Job(b, vv; path=path, freq=50)
        submit(job)
        return nothing
    end
    sim_osb(30, "$path_mpi")
    """
    run(`$(Peridynamics.MPI.mpiexec()) -n 2 $(Base.julia_cmd()) --project -e $(mpi_cmd)`)

    @test isdir(path_threads_vtk)
    @test isdir(path_mpi_vtk)
    vtk_files_threads = Peridynamics.find_vtk_files(path_threads_vtk)
    vtk_files_mpi = Peridynamics.find_vtk_files(path_mpi_vtk)
    @test length(vtk_files_mpi) == length(vtk_files_threads) == 3
    for i in eachindex(vtk_files_threads, vtk_files_mpi)
        res_threads = read_vtk(vtk_files_threads[i])
        res_mpi = read_vtk(vtk_files_mpi[i])
        for key in keys(res_threads)
            @test res_threads[key] ≈ res_mpi[key]
        end
    end

    rm(root; recursive=true, force=true)
end

@testitem "MPI-Threads comparison CCMaterial" begin
    root = joinpath(@__DIR__, "temp_mpi_threads_comparison_cc")
    path_threads = joinpath(root, "results_threads")
    path_threads_vtk = joinpath(path_threads, "vtk")
    path_mpi = joinpath(root, "results_mpi")
    path_mpi_vtk = joinpath(path_mpi, "vtk")

    function sim_cc(N::Int, path::String)
        l, Δx, δ, a = 1.0, 1/N, 3.015/N, 0.5
        pos, vol = uniform_box(l, l, 0.1l, Δx)
        ids = sortperm(pos[2,:])
        b = Body(CCMaterial(), pos[:, ids], vol[ids])
        material!(b; horizon=3.015Δx, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
        point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
        precrack!(b, :set_a, :set_b)
        point_set!(p -> p[2] > l/2-Δx, b, :set_top)
        point_set!(p -> p[2] < -l/2+Δx, b, :set_bottom)
        velocity_bc!(t -> -30, b, :set_bottom, :y)
        velocity_bc!(t -> 30, b, :set_top, :y)
        vv = VelocityVerlet(steps=100)
        job = Job(b, vv; path=path, freq=50)
        submit(job)
        return nothing
    end
    sim_cc(30, path_threads)

    mpi_cmd = """
    using Peridynamics
    function sim_cc(N::Int, path::String)
        l, Δx, δ, a = 1.0, 1/N, 3.015/N, 0.5
        pos, vol = uniform_box(l, l, 0.1l, Δx)
        ids = sortperm(pos[2,:])
        b = Body(CCMaterial(), pos[:, ids], vol[ids])
        material!(b; horizon=3.015Δx, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
        point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
        point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
        precrack!(b, :set_a, :set_b)
        point_set!(p -> p[2] > l/2-Δx, b, :set_top)
        point_set!(p -> p[2] < -l/2+Δx, b, :set_bottom)
        velocity_bc!(t -> -30, b, :set_bottom, :y)
        velocity_bc!(t -> 30, b, :set_top, :y)
        vv = VelocityVerlet(steps=100)
        job = Job(b, vv; path=path, freq=50)
        submit(job)
        return nothing
    end
    sim_cc(30, "$path_mpi")
    """
    run(`$(Peridynamics.MPI.mpiexec()) -n 2 $(Base.julia_cmd()) --project -e $(mpi_cmd)`)

    @test isdir(path_threads_vtk)
    @test isdir(path_mpi_vtk)
    vtk_files_threads = Peridynamics.find_vtk_files(path_threads_vtk)
    vtk_files_mpi = Peridynamics.find_vtk_files(path_mpi_vtk)
    @test length(vtk_files_mpi) == length(vtk_files_threads) == 3
    for i in eachindex(vtk_files_threads, vtk_files_mpi)
        res_threads = read_vtk(vtk_files_threads[i])
        res_mpi = read_vtk(vtk_files_mpi[i])
        for key in keys(res_threads)
            @test res_threads[key] ≈ res_mpi[key]
        end
    end

    rm(root; recursive=true, force=true)
end
