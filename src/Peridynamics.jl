module Peridynamics

using Base.Threads, Printf, LinearAlgebra, StaticArrays, PointNeighbors, ProgressMeter,
      WriteVTK, TimerOutputs, MPI, PrecompileTools
@static if Sys.islinux()
    using ThreadPinning
end

import LibGit2, Dates
import Polyester: @batch

# Material models
export BBMaterial, OSBMaterial, NOSBMaterial, CKIMaterial

# Systems related types
export NoCorrection, EnergySurfaceCorrection

# Discretization
export Body, point_set!, point_sets, failure_permit!, material!, velocity_bc!, velocity_ic!,
       forcedensity_bc!, precrack!, MultibodySetup, contact!, uniform_box, uniform_sphere,
       n_points

# Running simulations
export VelocityVerlet, DynamicRelaxation, Job, submit

# Pre processing
export read_inp

# Post processing and helpers
export read_vtk, process_each_export, mpi_isroot, force_mpi_run!, force_threads_run!,
       enable_mpi_timers!, disable_mpi_timers!, enable_mpi_progress_bars!,
       reset_mpi_progress_bars!, @mpitime, @mpiroot

function __init__()
    init_mpi()
    @static if Sys.islinux()
        mpi_run() || pinthreads(:cores; force=false)
    end
    BLAS.set_num_threads(1)
    return nothing
end

# """
#     AbstractMaterial

# Base for types which define materials.
# """
abstract type AbstractMaterial end

# """
#     AbstractSpatialSetup

# Base for types which define the spatial setup.
# """
abstract type AbstractSpatialSetup end

# """
#     AbstractBody{T<:AbstractMaterial} <: AbstractSpatialSetup

# Base for types which define a body.
# """
abstract type AbstractBody{T<:AbstractMaterial} <: AbstractSpatialSetup end

# """
#     AbstractMultibodySetup <: AbstractSpatialSetup

# Base for types which define a multibody setup.
# """
abstract type AbstractMultibodySetup <: AbstractSpatialSetup end

# """
#     AbstractParameterSetup

# Base for types which define a parameter setup.
# """
abstract type AbstractParameterSetup end

# """
#     AbstractPointParameters <: AbstractParameterSetup

# Base for types which define point parameters for specific materials.
# """
abstract type AbstractPointParameters <: AbstractParameterSetup end
abstract type AbstractParamSpec end

# """
#     AbstractTimeSolver

# Base for types which define a time solver.
# """
abstract type AbstractTimeSolver end

# """
#     AbstractJob

# Base for types which define a job.
# """
abstract type AbstractJob end
abstract type AbstractJobOptions end

# """
#     AbstractSystem

# Base for types which define a system.
# """
abstract type AbstractSystem end

# """
#     AbstractPredefinedCrack

# Base for types which define a predefined crack.
# """
abstract type AbstractPredefinedCrack end
abstract type AbstractBodyChunk{S<:AbstractSystem,T<:AbstractMaterial} end
abstract type AbstractParameterHandler <: AbstractParameterSetup end
abstract type AbstractChunkHandler end
abstract type AbstractDataHandler end
abstract type AbstractThreadsDataHandler <: AbstractDataHandler end
abstract type AbstractMPIDataHandler <: AbstractDataHandler end
abstract type AbstractThreadsBodyDataHandler{Sys,M,P,S} <: AbstractThreadsDataHandler end
abstract type AbstractThreadsMultibodyDataHandler <: AbstractThreadsDataHandler end
abstract type AbstractMPIBodyDataHandler{Sys,M,P,S} <: AbstractMPIDataHandler end
abstract type AbstractMPIMultibodyDataHandler <: AbstractMPIDataHandler end

# """
#     AbstractCorrection

# Base for types which define a correction algorithm suitable for materials that are a
# subtype of [`AbstractBondSystemMaterial`](@ref) and use the [`BondSystem`](@ref).
# """
abstract type AbstractCorrection end

# """
#     AbstractStorage

# Base for types which define a storage.
# """
abstract type AbstractStorage end

# """
#     AbstractStorage

# Base for types which define a condition.
# """
abstract type AbstractCondition end

# """
#     AbstractBondSystemMaterial{Correction} <: AbstractMaterial

# Base for types which define a material that uses the [`BondSystem`](@ref).
# """
abstract type AbstractBondSystemMaterial{Correction} <: AbstractMaterial end

# """
#     AbstractInteractionSystemMaterial <: AbstractMaterial

# Base for types which define a material that uses the [`InteractionSystem`](@ref).
# """
abstract type AbstractInteractionSystemMaterial <: AbstractMaterial end

include("auxiliary/function_arguments.jl")
include("auxiliary/io.jl")
include("auxiliary/logs.jl")
include("auxiliary/mpi.jl")

include("physics/boundary_conditions.jl")
include("physics/initial_conditions.jl")
include("physics/material_parameters.jl")
include("physics/fracture.jl")
include("physics/short_range_force_contact.jl")

include("discretization/point_generators.jl")
include("discretization/predefined_cracks.jl")
include("discretization/point_sets.jl")
include("discretization/body.jl")
include("discretization/multibody_setup.jl")
include("discretization/decomposition.jl")
include("discretization/chunk_handler.jl")
include("discretization/bond_system.jl")
include("discretization/bond_system_corrections.jl")
include("discretization/interaction_system.jl")
include("discretization/body_chunk.jl")

include("core/job.jl")
include("core/submit.jl")
include("core/parameter_handler.jl")
include("core/systems.jl")
include("core/materials.jl")
include("core/storages.jl")
include("core/time_solvers.jl")
include("core/halo_exchange.jl")
include("core/data_handler.jl")
include("core/threads_body_data_handler.jl")
include("core/threads_multibody_data_handler.jl")
include("core/mpi_body_data_handler.jl")
include("core/mpi_multibody_data_handler.jl")

include("time_solvers/velocity_verlet.jl")
include("time_solvers/dynamic_relaxation.jl")

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
