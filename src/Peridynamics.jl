module Peridynamics

using Printf
using ProgressMeter
using WriteVTK
using Base.Threads

export PointCloud, pcmerge
export PreCrack
export VelocityBC, VelocityIC, ForceDensityBC, PosDepVelBC
export VelocityVerlet, DynamicRelaxation, ExportSettings
export MultiMaterial, BondBasedMaterial, ContinuumBasedMaterial
export Contact, BodySetup
export PDSingleBodyAnalysis, PDContactAnalysis
export submit
export read_inp
export read_vtk

# Basics for peridynamics and contact problems
include("abstract_types.jl")
include("multi_material.jl")
include("spatial_discretization.jl")
# include("time_discretization.jl")
include("conditions.jl")
include("io.jl")
include("velocity_verlet.jl")
include("dynamic_relaxation.jl")
include("jobs.jl")
include("contact.jl")
include("utilities.jl")

# Material models
include("bond_based.jl")
include("continuum_based.jl")

# Helper modules
# Convert Abaqus mesh to PointCloud
include("AbaqusMeshConverter.jl")
using .AbaqusMeshConverter: read_inp

# Read VTK results
include("VtkReader.jl")
using .VtkReader: read_vtk

end
