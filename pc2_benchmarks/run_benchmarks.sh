#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=medium
#SBATCH --mem=0

module load lang       # loading the gateway module
module load JuliaHPC   # loading the latest JuliaHPC

julia setup.jl

respath="pc2_benchmarks"
rm -rf $respath
export_vtk=(true false)
threads=(1 4 8 16 32)
for t in ${threads[@]}; do
    for e in ${export_vtk[@]}; do
        julia --project -t $t bbvv.jl $respath "bbvv_t=${t}_e=$e" $e
    done
done

julia --project -t 8 create_benchmark_plots.jl $respath
