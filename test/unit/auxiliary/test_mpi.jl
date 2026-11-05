@testitem "force_mpi_run!" begin
    # initial setup of the MPI refs
    Peridynamics.MPI_RUN[] = false
    Peridynamics.MPI_RUN_FORCED[] = false

    Peridynamics.force_mpi_run!()
    @test Peridynamics.mpi_run() == true
    Peridynamics.init_mpi() # will do nothing because MPI_RUN was forced
    @test Peridynamics.mpi_run() == true

    # reset MPI refs
    Peridynamics.MPI_RUN[] = false
    Peridynamics.MPI_RUN_FORCED[] = false
end

@testitem "force_threads_run!" begin
    Peridynamics.force_threads_run!()
    @test Peridynamics.mpi_run() == false
    Peridynamics.init_mpi() # will do nothing because MPI_RUN was forced
    @test Peridynamics.mpi_run() == false
end

@testitem "@mpiroot error case" begin # lint-ok: the println is never evaluated
    @test_throws LoadError eval(:(@mpiroot :wrongkey println("Hi!")))
end

@testitem "mpi_barrier" begin
    # Test that mpi_barrier returns nothing
    @test mpi_barrier() === nothing
end
