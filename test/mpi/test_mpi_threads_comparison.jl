@testitem "MPI: results equal the threads results" tags=[:mpi, :simulation] setup=[Fixtures, MPIComparison] begin
    # The same simulations once with threads in this process and once with MPI in a subprocess,
    # compared file by file. The decomposition must not change the result. The cases are defined
    # once in `scripts/comparison_sims.jl` and run by `scripts/comparison.jl`.
    force_threads_run!()
    root = mktempdir()
    root_threads = joinpath(root, "threads")
    root_mpi = joinpath(root, "mpi")

    nranks = MPIComparison.COMPARISON_NRANKS
    for (name, _) in MPIComparison.COMPARISON_CASES
        path = MPIComparison.case_path(root_threads, name)
        MPIComparison.run_comparison_case(name, path; n_chunks=nranks)
    end

    @test Fixtures.run_mpi_script("comparison.jl", root_mpi; nranks)

    for (name, (_, n_files, tol)) in MPIComparison.COMPARISON_CASES
        @testset "$name" begin
            files_threads = Peridynamics.find_vtk_files(joinpath(MPIComparison.case_path(root_threads, name), "vtk"))
            files_mpi = Peridynamics.find_vtk_files(joinpath(MPIComparison.case_path(root_mpi, name), "vtk"))
            @test length(files_threads) == n_files
            @test length(files_mpi) == n_files
            for (ft, fm) in zip(files_threads, files_mpi)
                res_threads = read_vtk(ft)
                res_mpi = read_vtk(fm)
                @test keys(res_threads) == keys(res_mpi)
                for key in keys(res_threads)
                    if tol === nothing
                        @test res_threads[key] ≈ res_mpi[key]
                    else
                        @test maximum(abs.(res_threads[key] .- res_mpi[key])) < tol
                    end
                end
            end
        end
    end
end
