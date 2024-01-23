using Test
using SafeTestsets

@testset "Peridynamics" verbose=true begin
    # unit tests
    @safetestset "point cloud" begin include(joinpath("unit","discretizations","test_point_cloud.jl")) end
    @safetestset "decomposition" begin include(joinpath("unit","discretizations","test_decomposition.jl")) end
    @safetestset "point bond discretization" begin include(joinpath("unit","discretizations","test_point_bond_discretization.jl")) end
    @safetestset "chunk local handler" begin include(joinpath("unit","threads_core","test_chunk_local_handler.jl")) end
end
