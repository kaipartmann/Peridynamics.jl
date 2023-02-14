module Peridynamics

using Printf
using ProgressMeter
using WriteVTK
using Base.Threads

export PointCloud
export PreCrack
export VelocityBC, VelocityIC, ForceDensityBC, PosDepVelBC
export TimeDiscretization, calc_stable_user_timestep, ExportSettings
export BondBasedMaterial
export Contact, BodySetup
export PDSingleBodyAnalysis, PDContactAnalysis
export submit
export read_inp
export read_vtk

# Basics for peridynamics and contact problems
include("abstract_types.jl")
include("spatial_discretization.jl")
include("time_discretization.jl")
include("conditions.jl")
include("io.jl")
include("jobs.jl")
include("contact.jl")
include("utilities.jl")

# Material models
include("bond_based.jl")

# Helper modules
# Convert Abaqus mesh to PointCloud
include("AbaqusMeshConverter.jl")
using .AbaqusMeshConverter: read_inp

# Read VTK results
include("VtkReader.jl")
using .VtkReader: read_vtk

end
