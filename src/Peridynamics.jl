module Peridynamics

using Base.Threads, Printf, LinearAlgebra, StaticArrays, ProgressMeter, WriteVTK,
      TimerOutputs, MPI

export PointCloud, BBMaterial, Body, point_set!, material!, velocity_bc!, velocity_ic!,
       forcedensity_bc!, forcedensity_ic!

const MPI_INITIALIZED = Ref(false)
const MPI_RANK = Ref(-1)
const MPI_SIZE = Ref(-1)
const MPI_SIM = Ref(true)
const TO = TimerOutput()
@inline mpi_comm() = MPI.COMM_WORLD
@inline mpi_rank() = MPI_RANK[]
@inline mpi_nranks() = MPI_SIZE[]
@inline mpi_sim() = MPI_SIM[]

const FIND_POINTS_ALLOWED_SYMBOLS = (:x, :y, :z, :p)
const SYMBOL_TO_DIM = Dict(:x => 0x1, :y => 0x2, :z => 0x3)
const DimensionSpec = Union{Integer,Symbol}

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

abstract type AbstractDiscretization end

abstract type AbstractPredefinedCrack end

abstract type AbstractDataHandler end

abstract type AbstractStorage end

abstract type AbstractCondition end

include("conditions/boundary_conditions.jl")
include("conditions/initial_conditions.jl")
include("conditions/condition_checks.jl")

include("discretizations/predefined_cracks.jl")
include("discretizations/point_sets.jl")
include("discretizations/find_points.jl")
include("discretizations/body.jl")

include("discretizations/point_cloud.jl")
include("discretizations/decomposition.jl")
include("discretizations/point_bond_discretization.jl")

include("materials/bond_based.jl")

include("threads_core/chunk_local_handler.jl")
include("threads_core/threads_halo_exchange.jl")
include("threads_core/threads_data_handler.jl")

include("auxiliary/function_arguments.jl")

end
