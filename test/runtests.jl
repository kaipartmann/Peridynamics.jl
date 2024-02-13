using Test
using SafeTestsets

@testset "Peridynamics" verbose=true begin
    # unit tests
    @safetestset "body" begin include("unit/discretizations/test_body.jl") end
    @safetestset "find points" begin include("unit/discretizations/test_find_points.jl") end
    @safetestset "decomposition" begin include("unit/discretizations/test_decomposition.jl") end
    @safetestset "function arguments" begin include("unit/auxiliary/test_function_arguments.jl") end
    @safetestset "velocity verlet" begin include("unit/time_solvers/test_velocity_verlet.jl") end
end
