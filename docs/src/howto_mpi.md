# Simulations with MPI

The package is designed so that the same core functions are used and only a small backend handles the differences between MPI or multithreading.
This means, the development goal was:

**Code that runs with multithreading should also work with MPI without changes!**

However, currently not all features are supported with MPI. A table with an overview is shown below.

#### Currently supported features:

| Job type | MPI |
|:---|:---|
| `Job(::Body, ::VelocityVerlet)` | ✅ | 
| `Job(::Body, ::DynamicRelaxation)` | ✅ | 
| `Job(::MultibodySetup, ::VelocityVerlet)` | ❌ | 

## Setting up simulations for MPI

If a script containing a simulation runs with multithreading and the features are supported with MPI, then this same script can be run with:
```bash
mpiexecjl -n <number of ranks> julia --project path/to/script.jl
```
Please refer to the [`MPI.jl` documentation of `mpiexecjl`](https://juliaparallel.org/MPI.jl/latest/usage/#Installation) for installation and setup instructions.

Furthermore, there are helper functions that improve the setup of MPI simulations, such as [`enable_mpi_timers!`](@ref), [`@mpiroot`](@ref), [`@mpitime`](@ref), or [`mpi_isroot`](@ref).

