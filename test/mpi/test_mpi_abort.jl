@testitem "MPI: submit aborts on a non-synchronized error" tags=[:mpi] setup=[Fixtures] begin
    # An error that is not a `NaNError` is not synchronized across the ranks, so `submit` writes
    # it to the logfile and aborts MPI instead of risking a deadlock. The process must exit
    # non-zero. This cannot share a process with any other MPI check.
    path = mktempdir()
    # somehow on Julia ≤ 1.10 the abort does not result in a non-zero exit status
    if VERSION > v"1.10"
        @test Fixtures.run_mpi_script("abort.jl", path; nranks=2, expect_success=false)
        @test occursin("some weird error occurred!", read(joinpath(path, "logfile.log"), String))
    end
end
