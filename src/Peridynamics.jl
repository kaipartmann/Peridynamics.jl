module Peridynamics

using FileIO
using Printf
using ProgressMeter
using WriteVTK
using StaticArrays
using LinearAlgebra
using AbaqusReader: abaqus_read_mesh
using Base.Threads

export PointCloud, read_inp
export PreCrack
export VelocityBC, VelocityIC, ForceDensityBC, PosDepVelBC
export TimeDiscretization, calc_stable_user_timestep, ExportSettings
export BondBasedMaterial
export Contact, BodySetup
export PDSingleBodyAnalysis, PDContactAnalysis
export submit

# Base modules for peridynamics and contact problems
include("PeridynamicsBase.jl")
include("ContactBase.jl")

# Convert Abaqus mesh to PointCloud
include("AbaqusMeshConverter.jl")

# Material models
include("BondBased.jl")

end
