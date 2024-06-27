@testitem "force_mpi_run!" begin
    Peridynamics.force_mpi_run!()
    @test Peridynamics.mpi_run() == true
    Peridynamics.set_mpi_run!(false) # will do nothing because MPI_RUN was forced
    @test Peridynamics.mpi_run() == true
end

@testitem "force_threads_run!" begin
    Peridynamics.force_threads_run!()
    @test Peridynamics.mpi_run() == false
    Peridynamics.set_mpi_run!(true) # will do nothing because MPI_RUN was forced
    @test Peridynamics.mpi_run() == false
end
