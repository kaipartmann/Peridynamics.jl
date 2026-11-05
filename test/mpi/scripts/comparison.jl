# MPI side of the MPI-vs-threads comparison. Run with
#     mpiexec -n 2 julia --project=<Peridynamics> test/mpi/scripts/comparison.jl <root>
# and find the results of every case of `comparison_sims.jl` below `<root>`.
using Peridynamics
include(joinpath(@__DIR__, "comparison_sims.jl"))
const ROOT = ARGS[1]
for (name, _) in COMPARISON_CASES
    run_comparison_case(name, case_path(ROOT, name))
end
