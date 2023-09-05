#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=medium
#SBATCH --mem=0

module load lang       # loading the gateway module
module load JuliaHPC   # loading the latest JuliaHPC

julia setup.jl

respath="pc2_bench_v1"
rm -rf $respath
export_vtk=(true)
threads=(1 8 16 32)
for t in ${threads[@]}; do
    for e in ${export_vtk[@]}; do
        julia --project -t $t bbvv.jl $respath "bbvv_t=${t}_e=$e" $e
    done
done
