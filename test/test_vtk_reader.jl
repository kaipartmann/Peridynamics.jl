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
        vtk["Damage", VTKPointData()] = damage
        vtk["Displacement", VTKPointData()] = displacement
        vtk["Time", VTKFieldData()] = time
    end
    result = read_vtk(name)

    @test typeof(result) == Dict{String, VecOrMat{Float64}}
    @test result["Position"] == position
    @test result["Time"] == [time]
    @test result["Damage"] == damage
    @test result["Displacement"] == displacement

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
    random_point_data = rand(n_points)
    random_field_data =  [1e-5, 1.0, 5]
    vtk_grid(bname, position, cells) do vtk
        vtk["Time", VTKFieldData()] = time
        vtk["Damage", VTKPointData()] = damage
        vtk["RandomPointData", VTKPointData()] = random_point_data
        vtk["RandomFieldData", VTKFieldData()] = random_field_data
    end
    result = read_vtk(name)

    @test typeof(result) == Dict{String, VecOrMat{Float64}}
    @test result["Position"] == position
    @test result["Time"] == [time]
    @test result["Damage"] == damage
    @test result["RandomPointData"] == random_point_data
    @test result["RandomFieldData"] == random_field_data

    rm(name)
end

@testset "wrong file type" begin
    @test_throws AssertionError read_vtk("something.wrong")
end

@testset "corrupt file <AppendedData encoding=\"raw\">" begin
    bname = joinpath(@__DIR__, "vtk_test_2")
    name = bname * ".vtu"
    isfile(name) && rm(name)
    n_points = 10
    position = rand(3, n_points)
    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (j,)) for j in 1:n_points]
    damage = rand(n_points)
    time = rand() * 1000
    random_point_data = rand(n_points)
    random_field_data =  [1e-5, 1.0, 5]
    vtk_grid(bname, position, cells) do vtk
        vtk["Time", VTKFieldData()] = time
        vtk["Damage", VTKPointData()] = damage
        vtk["RandomPointData", VTKPointData()] = random_point_data
        vtk["RandomFieldData", VTKFieldData()] = random_field_data
    end
    file_raw = read(name, String)
    file_raw_new = replace(file_raw,"<AppendedData encoding=\"raw\">" => "<>",)
    open(name, "w") do io
        write(io, file_raw_new)
    end

    @test_throws ErrorException read_vtk(name)

    rm(name)
end

@testset "corrupt file offset marker" begin
    bname = joinpath(@__DIR__, "vtk_test_2")
    name = bname * ".vtu"
    isfile(name) && rm(name)
    n_points = 10
    position = rand(3, n_points)
    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (j,)) for j in 1:n_points]
    damage = rand(n_points)
    time = rand() * 1000
    random_point_data = rand(n_points)
    random_field_data =  [1e-5, 1.0, 5]
    vtk_grid(bname, position, cells) do vtk
        vtk["Time", VTKFieldData()] = time
        vtk["Damage", VTKPointData()] = damage
        vtk["RandomPointData", VTKPointData()] = random_point_data
        vtk["RandomFieldData", VTKFieldData()] = random_field_data
    end
    file_raw = read(name, String)
    file_raw_new = replace(file_raw, "_" => "")
    open(name, "w") do io
        write(io, file_raw_new)
    end

    @test_throws ErrorException read_vtk(name)

    rm(name)
end

@testset "corrupt file </AppendedData>" begin
    bname = joinpath(@__DIR__, "vtk_test_2")
    name = bname * ".vtu"
    isfile(name) && rm(name)
    n_points = 10
    position = rand(3, n_points)
    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (j,)) for j in 1:n_points]
    damage = rand(n_points)
    time = rand() * 1000
    random_point_data = rand(n_points)
    random_field_data =  [1e-5, 1.0, 5]
    vtk_grid(bname, position, cells) do vtk
        vtk["Time", VTKFieldData()] = time
        vtk["Damage", VTKPointData()] = damage
        vtk["RandomPointData", VTKPointData()] = random_point_data
        vtk["RandomFieldData", VTKFieldData()] = random_field_data
    end
    file_raw = read(name, String)
    file_raw_new = replace(file_raw, "\n  </AppendedData>" => "<AppendedData>")
    open(name, "w") do io
        write(io, file_raw_new)
    end

    @test_throws ErrorException read_vtk(name)

    rm(name)
end
