import Pkg; Pkg.activate(@__DIR__)

using BenchmarkTools
using Peridynamics

const SUITE = BenchmarkGroup()

include(joinpath(@__DIR__, "mode_i.jl"))

SUITE["mode_i"] = BenchmarkGroup()
SUITE["mode_i"]["BBMaterial, 40"] = @benchmarkable submit(mode_i(BBMaterial(), 40))
SUITE["mode_i"]["OSBMaterial, 40"] = @benchmarkable submit(mode_i(OSBMaterial(), 40))
SUITE["mode_i"]["CCMaterial, 40"] = @benchmarkable submit(mode_i(CCMaterial(), 40))
