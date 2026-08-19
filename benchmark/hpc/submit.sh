#!/bin/bash
# Slurm template for the MPI scaling benchmark, submitted from the root of the checkout with
#
#     sbatch benchmark/hpc/submit.sh
#
# The SBATCH header has to be adjusted once per cluster.
#SBATCH --job-name=pd-scaling
#SBATCH --time=08:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kai.partmann@uni-siegen.de

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PROJECT="$REPO/benchmark/hpc"
LABEL="${SLURM_CLUSTER_NAME:-$(hostname)}"

#-- environment, once and serially ---------------------------------------------------------
# `scaling.jl` does no package management of its own: a few hundred ranks resolving the same
# environment simultaneously is how a shared depot gets corrupted.
julia --project="$PROJECT" -e "
    import Pkg
    Pkg.develop(Pkg.PackageSpec(path=\"$REPO\"))
    Pkg.instantiate()
    using Peridynamics, MPI
"

#-- smoke test -----------------------------------------------------------------------------
# Seconds, and it fails immediately if the launcher, the MPI build or the checkout is wrong.
mpiexecjl -n 2 julia --project="$PROJECT" "$PROJECT/scaling.jl" --size=smoke --label="$LABEL"

#-- the sweep ------------------------------------------------------------------------------
# One invocation per point of the scaling curve. Every rank count has to fit into the
# allocation above, and the first one is the reference the speedups are measured against, so it
# has to be a configuration the problem actually fits into.
for n in 8 16 32 64 128 256; do
    mpiexecjl -n "$n" julia --project="$PROJECT" "$PROJECT/scaling.jl" \
        --size=large --mat=BB --label="$LABEL"

    # If `mpiexecjl` writes its result and then does not return, that is the hydra launcher and
    # not the measurement. Use Slurm's own launcher instead, which tears its tasks down:
    #
    #     srun --ntasks="$n" julia --project="$PROJECT" "$PROJECT/scaling.jl" \
    #         --size=large --mat=BB --label="$LABEL"
done

#-- afterwards -----------------------------------------------------------------------------
# Copy `benchmark/hpc/results/` back and run
#     julia --project=benchmark/hpc benchmark/hpc/report.jl
