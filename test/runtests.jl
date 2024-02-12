using Test
using SafeTestsets

@testset "Peridynamics" verbose=true begin
    # unit tests
    @safetestset "point cloud" begin include("unit/discretizations/test_point_cloud.jl") end
    @safetestset "decomposition" begin include("unit/discretizations/test_decomposition.jl") end
    @safetestset "point bond discretization" begin include("unit/discretizations/test_point_bond_discretization.jl") end
    @safetestset "chunk local handler" begin include("unit/threads_core/test_chunk_local_handler.jl") end
    @safetestset "function arguments" begin include("unit/auxiliary/test_function_arguments.jl") end
    @safetestset "body" begin include("unit/discretizations/test_body.jl") end
    @safetestset "find points" begin include("unit/discretizations/test_find_points.jl") end
    @safetestset "velocity verlet" begin include("unit/time_solvers/test_velocity_verlet.jl") end
end
