using Peridynamics
using Test
using SafeTestsets

@testset "Peridynamics.jl" begin
    @safetestset "Types" begin include("types.jl") end
    @safetestset "Bonds" begin include("bonds.jl") end
    @safetestset "BondBased" begin include("bondbased.jl") end
    @safetestset "Symmetry" begin include("symmetry.jl") end
    @safetestset "AbaqusMeshConverter" begin include("meshconverter.jl") end
end
