module Peridynamics

using Base.Threads, Printf, LinearAlgebra, StaticArrays, ProgressMeter, WriteVTK,
      TimerOutputs, MPI

export BBMaterial, CKIMaterial, Body, point_set!, failure_permit!, material!, velocity_bc!,
       velocity_ic!, forcedensity_bc!, precrack!, VelocityVerlet, MultibodySetup, contact!,
       Job, read_vtk, uniform_box, submit

const MPI_INITIALIZED = Ref(false)
const MPI_RANK = Ref(-1)
const MPI_SIZE = Ref(-1)
const MPI_SIM = Ref(true)
const QUIET = Ref(false)
@inline mpi_comm() = MPI.COMM_WORLD
@inline mpi_rank() = MPI_RANK[]
@inline mpi_nranks() = MPI_SIZE[]
@inline mpi_sim() = MPI_SIM[]
@inline quiet() = QUIET[]
@inline function set_quiet(b::Bool)
    QUIET[] = b
    return nothing
end

const TO = TimerOutput()

const FIND_POINTS_ALLOWED_SYMBOLS = (:x, :y, :z, :p)
const SYMBOL_TO_DIM = Dict(:x => 0x1, :y => 0x2, :z => 0x3)
const ELASTIC_KWARGS = (:E, :nu)
const FRAC_KWARGS = (:Gc, :epsilon_c)
const DEFAULT_POINT_KWARGS = (:horizon, :rho, ELASTIC_KWARGS..., FRAC_KWARGS...)
const CONTACT_KWARGS = (:radius, :sc)
const EXPORT_KWARGS = (:path, :freq, :write)
const JOB_KWARGS = (EXPORT_KWARGS...,)
const SUBMIT_KWARGS = (:quiet,)

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
abstract type AbstractPointParameters end

abstract type AbstractTimeSolver end

abstract type AbstractDiscretization end

abstract type AbstractPredefinedCrack end

abstract type AbstractDataHandler end

abstract type AbstractStorage end

abstract type AbstractCondition end

include("conditions/boundary_conditions.jl")
include("conditions/initial_conditions.jl")
include("conditions/condition_checks.jl")

include("discretizations/point_generators.jl")
include("discretizations/predefined_cracks.jl")
include("discretizations/find_points.jl")
include("discretizations/body.jl")
include("discretizations/contact.jl")
include("discretizations/multibody_setup.jl")
const SpatialSetup = Union{Body,MultibodySetup}
include("discretizations/decomposition.jl")
include("discretizations/bond_discretization.jl")
include("discretizations/chunk_handler.jl")
include("discretizations/body_chunk.jl")

include("time_solvers/velocity_verlet.jl")

include("physics/material_interface.jl")
include("physics/material_parameters.jl")
include("physics/fracture.jl")
include("physics/bond_based.jl")
include("physics/continuum_kinematics_inspired.jl")

include("auxiliary/function_arguments.jl")
include("auxiliary/logs.jl")
include("auxiliary/io.jl")

include("core/job.jl")
include("core/submit.jl")
include("core/submit_threads.jl")
include("core/submit_mpi.jl")

include("VtkReader/VtkReader.jl")
using .VtkReader

include("AbaqusMeshConverter/AbaqusMeshConverter.jl")
using .AbaqusMeshConverter

end
