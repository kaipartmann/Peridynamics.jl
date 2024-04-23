module Peridynamics

using Base.Threads, Printf, LinearAlgebra, StaticArrays, NearestNeighbors, ProgressMeter,
      WriteVTK, TimerOutputs, MPI, PrecompileTools
@static if Sys.islinux()
    using ThreadPinning
end

import LibGit2, Dates

export BBMaterial, CKIMaterial, NOSBMaterial, OSBMaterial, Body, point_set!,
       failure_permit!, material!, velocity_bc!, velocity_ic!, forcedensity_bc!, precrack!,
       VelocityVerlet, MultibodySetup, contact!, Job, read_vtk, uniform_box, submit,
       process_each_export, mpi_isroot, force_mpi_run!, force_threads_run!,
       enable_mpi_timers!, disable_mpi_timers!, @mpitime, @rootdo

function __init__()
    init_mpi()
    @static if Sys.islinux()
        mpi_run() || pinthreads(:cores; force=false)
    end
    BLAS.set_num_threads(1)
    return nothing
end

abstract type AbstractMaterial end
abstract type AbstractSpatialSetup end
abstract type AbstractBody{T<:AbstractMaterial} <: AbstractSpatialSetup end
abstract type AbstractMultibodySetup{T<:AbstractMaterial} <: AbstractSpatialSetup end
abstract type AbstractPointParameters end
abstract type AbstractTimeSolver end
abstract type AbstractJob end
abstract type AbstractOptions end
abstract type AbstractSystem end
abstract type AbstractPredefinedCrack end
abstract type AbstractBodyChunk{System,Material} end
abstract type AbstractParameterHandler{N} end
abstract type AbstractChunkHandler end
abstract type AbstractDataHandler end
abstract type AbstractThreadsDataHandler <: AbstractDataHandler end
abstract type AbstractMPIDataHandler <: AbstractDataHandler end
abstract type AbstractCorrectionHandler end
abstract type AbstractStorage end
abstract type AbstractCondition end

include("auxiliary/function_arguments.jl")
include("auxiliary/io.jl")
include("auxiliary/logs.jl")
include("auxiliary/mpi.jl")

include("conditions/boundary_conditions.jl")
include("conditions/initial_conditions.jl")
include("conditions/condition_checks.jl")

include("discretization/point_generators.jl")
include("discretization/predefined_cracks.jl")
include("discretization/point_sets.jl")
include("discretization/body.jl")
include("discretization/contact.jl")
include("discretization/multibody_setup.jl")
include("discretization/decomposition.jl")
include("discretization/chunk_handler.jl")
include("discretization/bond_system.jl")
include("discretization/body_chunk.jl")

include("core/job.jl")
include("core/submit.jl")
include("core/parameter_handler.jl")
include("core/systems.jl")
include("core/materials.jl")
include("core/storages.jl")
include("core/time_solvers.jl")
include("core/halo_exchange.jl")
include("core/threads_data_handler.jl")
include("core/mpi_data_handler.jl")

include("time_solvers/velocity_verlet.jl")

include("physics/force_density.jl")
include("physics/material_parameters.jl")
include("physics/fracture.jl")
include("physics/bond_based.jl")
include("physics/continuum_kinematics_inspired.jl")
include("physics/ordinary_state_based.jl")
include("physics/correspondence.jl")

include("VtkReader/VtkReader.jl")
using .VtkReader

include("AbaqusMeshConverter/AbaqusMeshConverter.jl")
using .AbaqusMeshConverter

include("auxiliary/process_each_export.jl")

try
    include("auxiliary/precompile_workload.jl")
catch err
    @error "precompilation errored\n" exception=err
end

end
