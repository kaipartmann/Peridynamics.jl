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
export read_vtk, SimResult

# Basics for peridynamics and contact problems
include("peridynamics_base.jl")
include("contact_base.jl")

# Material models
include("bond_based.jl")

# Convert Abaqus mesh to PointCloud
include("AbaqusMeshConverter.jl")
using .AbaqusMeshConverter: read_inp

# Read VTK results
include("VtkReader.jl")
using .VtkReader: read_vtk, SimResult

end
