module Peridynamics

using Base.Threads, Printf, LinearAlgebra, StaticArrays, ProgressMeter, WriteVTK,
    TimerOutputs, MPI

export PointCloud, BBMaterial

const MPI_INITIALIZED = Ref(false)
const MPI_RANK = Ref(-1)
const MPI_SIZE = Ref(-1)
const MPI_SIM = Ref(true)
const TO = TimerOutput()
@inline mpi_comm() = MPI.COMM_WORLD
@inline mpi_rank() = MPI_RANK[]
@inline mpi_nranks() = MPI_SIZE[]
@inline mpi_sim() = MPI_SIM[]

function __init__()
    if MPI_INITIALIZED[]
        return nothing
    end
    MPI.Init(finalize_atexit=true)
    MPI_RANK[] = MPI.Comm_rank(MPI.COMM_WORLD)
    MPI_SIZE[] = MPI.Comm_size(MPI.COMM_WORLD)
    MPI_INITIALIZED[] = true
    if !haskey(ENV, "MPI_LOCALRANKID")
        MPI_SIM[] = false
    end
    return nothing
end

include("basic_types.jl")

include(joinpath("discretizations", "point_cloud.jl"))
include(joinpath("discretizations", "decomposition.jl"))
include(joinpath("discretizations", "point_bond_discretization.jl"))

include(joinpath("materials", "bond_based.jl"))

include(joinpath("threads_core", "chunk_local_handler.jl"))
include(joinpath("threads_core", "threads_data_handler.jl"))

end
