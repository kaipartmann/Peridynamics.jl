using Peridynamics
using Test
using SafeTestsets

@testset "Peridynamics.jl" verbose=true begin
    @safetestset "Types" begin include("types.jl") end
    @safetestset "Bonds" begin include("bonds.jl") end
    @safetestset "IO" begin include("io.jl") end
    @safetestset "BondBased" begin include("bondbased.jl") end
    @safetestset "Symmetry" begin include("symmetry.jl") end
    @safetestset "AbaqusMeshConverter" begin include("meshconverter.jl") end
    @safetestset "Contact" begin include("contact.jl") end
    @safetestset "Utilities" begin include("utilities.jl") end
    @safetestset "TimeIntegration" begin include("timeintegration.jl") end
end
