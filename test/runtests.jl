using Peridynamics
using Test
using SafeTestsets

@testset "Peridynamics.jl" begin
    @safetestset "Bonds" begin include("bonds.jl") end
end
