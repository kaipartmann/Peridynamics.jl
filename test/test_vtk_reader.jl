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

@testset "show SimResult" begin
    position = rand(3, 5)
    time = 1.0
    damage = rand(5)
    displacement = position .+ 0.1
    sr = SimResult(position, time, damage, displacement)
    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (20, 40)), "text/plain", sr)
    msg_sr = String(take!(io))
    @test msg_sr == string(
        "SimResult with fields:\n",
        " position:     3×5 Matrix{Float64}\n",
        " time:         Float64\n",
        " damage:       5-element Vector{Float64}\n",
        " displacement: 3×5 Matrix{Float64}",
    )
end
