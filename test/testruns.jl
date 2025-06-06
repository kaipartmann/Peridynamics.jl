@testitem "MPI.jl MWE" begin
    mpi_cmd = """
    @show Base.active_project()
    using Peridynamics
    sleep(0.1 * Peridynamics.mpi_rank())
    println("Hello from process ", Peridynamics.mpi_rank())
    """
    run(`$(Peridynamics.MPI.mpiexec()) -n 2 $(Base.julia_cmd()) --project -e $(mpi_cmd)`)
end


@testitem "MWE loading Peridynamics" begin
    mpi_cmd = """
    @show Base.active_project()
    using Peridynamics
    """
    run(`$(Base.julia_cmd()) --project -e $(mpi_cmd)`)
end
