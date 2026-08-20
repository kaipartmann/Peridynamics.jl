# Threads side of the MPI-vs-threads comparison: the same simulations the `mpiexec` subprocess
# runs, see `scripts/comparison_sims.jl`.
@testmodule MPIComparison begin
    using Peridynamics
    include(joinpath(@__DIR__, "scripts", "comparison_sims.jl"))
end
