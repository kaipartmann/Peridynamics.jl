module Peridynamics

using Base.Threads, Printf, LinearAlgebra, StaticArrays, ProgressMeter, WriteVTK,
    TimerOutputs, MPI

export PointCloud, BBMaterial, Body, point_set!, material!

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

abstract type AbstractJob end

abstract type AbstractMaterial end
abstract type AbstractMaterialConfig end
abstract type MaterialHandler <: AbstractMaterialConfig end
abstract type Material <: AbstractMaterialConfig end

abstract type AbstractTimeSolver end

# abstract type AbstractDiscretizationConfig end
abstract type AbstractDiscretization end

abstract type AbstractPredefinedCrack end

abstract type AbstractDataHandler end
# abstract type MPIDataHandler <: AbstractDataHandler end
# abstract type ThreadsDataHandler <: AbstractDataHandler end

abstract type AbstractStorage end

abstract type AbstractCondition end
abstract type AbstractBoundaryCondition <: AbstractCondition end
abstract type AbstractSingleValueBC <: AbstractBoundaryCondition end

abstract type AbstractInitialCondition <: AbstractCondition end
abstract type AbstractSingleValueIC <: AbstractInitialCondition end

abstract type AbstractPointSetHandler end


include("discretizations/point_sets.jl")
include("discretizations/body.jl")
include("discretizations/find_points.jl")

include("discretizations/point_cloud.jl")
include("discretizations/decomposition.jl")
include("discretizations/point_bond_discretization.jl")

include("conditions/boundary_conditions.jl")
include("conditions/initial_conditions.jl")

include("materials/bond_based.jl")

include("threads_core/chunk_local_handler.jl")
include("threads_core/threads_halo_exchange.jl")
include("threads_core/threads_data_handler.jl")

include("auxiliary/function_arguments.jl")

end
