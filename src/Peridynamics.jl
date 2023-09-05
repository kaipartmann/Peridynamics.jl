module Peridynamics

using Printf
using ProgressMeter
using WriteVTK
using Base.Threads
using TimerOutputs
using ThreadPinning

export PointCloud, PreCrack, pcmerge
export VelocityBC, VelocityIC, ForceDensityBC, PosDepVelBC
# export VelocityVerlet, DynamicRelaxation, ExportSettings
export VelocityVerlet, ExportSettings
# export MultiMaterial, BondBasedMaterial, ContinuumBasedMaterial
export MultiMaterial, BondBasedMaterial
# export Contact, BodySetup
# export PDSingleBodyAnalysis, PDContactAnalysis, submit
export PDSingleBodyAnalysis, submit
export read_inp, read_vtk

const TO = TimerOutput()

# default ThreadPinning setting
pinthreads(:cores; force=false)

# Basics for peridynamics and contact problems
include("abstract_types.jl")
include("multi_material.jl")
include("spatial_discretization.jl")
include("conditions.jl")
include("io.jl")
include("pdproblem.jl")
include("velocity_verlet.jl")
# include("dynamic_relaxation.jl")
include("jobs.jl")
# include("contact.jl")
include("utilities.jl")

# Material models
include("material_interface.jl")
include("bond_based.jl")
# include("continuum_based.jl")

# Helper modules
# Convert Abaqus mesh to PointCloud
include("AbaqusMeshConverter.jl")
using .AbaqusMeshConverter: read_inp

# Read VTK results
include("VtkReader.jl")
using .VtkReader: read_vtk

end
