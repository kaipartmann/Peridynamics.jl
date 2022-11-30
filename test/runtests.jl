using Peridynamics
using Test
using SafeTestsets

@testset "Peridynamics.jl" verbose=true begin
    @safetestset "spatial_discretization" begin include("test_spatial_discretization.jl") end
    @safetestset "time_discretization" begin include("test_time_discretization.jl") end
    @safetestset "conditions" begin include("test_conditions.jl") end
    @safetestset "io" begin include("test_io.jl") end
    @safetestset "jobs" begin include("test_jobs.jl") end
    @safetestset "contact" begin include("test_contact.jl") end
    @safetestset "utilities" begin include("test_utilities.jl") end
    @safetestset "bond_based" begin include("test_bond_based.jl") end
    @safetestset "AbaqusMeshConverter" begin include("test_abaqus_mesh_converter.jl") end
    @safetestset "VtkReader" begin include("test_vtk_reader.jl") end
    @safetestset "MWE symmetry" begin include("test_mwe_symmetry.jl") end
end
