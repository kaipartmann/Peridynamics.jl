using Peridynamics
using Test

# @test Peridynamics.switch_to_threads() == true
Peridynamics.mpi_rank() == 0 ? println(Peridynamics.switch_to_threads()) : nothing
