using Test
using WriteVTK
using Peridynamics

@testset "read complete results" begin
    bname = joinpath(@__DIR__, "vtk_test_1")
    name = bname * ".vtu"
    isfile(name) && rm(name)
    n_points = 10
    position = rand(3, n_points)
    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (j,)) for j in 1:n_points]
    damage = rand(n_points)
    displacement = rand(3, n_points)
    time = rand() * 1000
    vtk_grid(bname, position, cells) do vtk
        vtk["damage", VTKPointData()] = damage
        vtk["displacement", VTKPointData()] = displacement
        vtk["time", VTKFieldData()] = time
    end
    result = read_vtk(name)

    @test typeof(result) == SimResult
    @test result.position == position
    @test result.time == time
    @test result.damage == damage
    @test result.displacement == displacement

    rm(name)
end

@testset "read incomplete results" begin
    bname = joinpath(@__DIR__, "vtk_test_2")
    name = bname * ".vtu"
    isfile(name) && rm(name)
    n_points = 10
    position = rand(3, n_points)
    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (j,)) for j in 1:n_points]
    damage = rand(n_points)
    time = rand() * 1000
    vtk_grid(bname, position, cells) do vtk
        vtk["damage", VTKPointData()] = damage
        vtk["random", VTKPointData()] = rand(n_points)
        vtk["time", VTKFieldData()] = time
    end
    result = read_vtk(name)

    @test typeof(result) == SimResult
    @test result.position == position
    @test result.time == time
    @test result.damage == damage
    @test result.displacement == Array{Float64, 2}(undef, 0, 0)

    rm(name)
end

@testset "wrong file type" begin
    @test_throws AssertionError read_vtk("something.wrong")
end
