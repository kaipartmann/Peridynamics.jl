using Peridynamics
using Test
using SafeTestsets

@testset "Peridynamics.jl" begin
    @safetestset "Types" begin include("types.jl") end
    @safetestset "Bonds" begin include("bonds.jl") end
    @safetestset "Symmetry" begin include("symmetry.jl") end
end
