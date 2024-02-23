

@testitem "read complete results" begin
    using Peridynamics.WriteVTK

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

#-- read incomplete results
@testitem "read incomplete results" begin
    using Peridynamics.WriteVTK

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

@testitem "wrong file type" begin
    @test_throws AssertionError read_vtk("something.wrong")
end

@testitem "corrupt file raw encoding" begin
    using Peridynamics.WriteVTK

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

@testitem "corrupt file offset marker" begin
    using Peridynamics.WriteVTK

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

#-- corrupt file </AppendedData>
@testitem "corrupt file appended data" begin
    using Peridynamics.WriteVTK

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
