#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=medium
#SBATCH --mem=0

julia setup.jl

respath="benchmark_results"
export_vtk=(false true)
threads=(1 2 4 8 16 32 64)
for t in ${threads[@]}; do
    for e in ${export_vtk[@]}; do
        julia --project -t $t bbvv.jl $respath "bbvv_t=${t}_e=$e" $e
        julia --project -t $t bbdr.jl $respath "bbdr_t=${t}_e=$e" $e
        julia --project -t $t cpdvv.jl $respath "cpdvv_t=${t}_e=$e" $e
        julia --project -t $t bbmmvv.jl $respath "bbmmvv_t=${t}_e=$e" $e
        julia --project -t $t contact.jl $respath "contact_t=${t}_e=$e" $e
    done
done
