

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
    integer_test = rand(1:100, n_points)
    time = rand() * 1000
    vtk_grid(bname, position, cells) do vtk
        vtk["Damage", VTKPointData()] = damage
        vtk["Displacement", VTKPointData()] = displacement
        vtk["integer_test", VTKPointData()] = integer_test
        vtk["Time", VTKFieldData()] = time
    end
    result = read_vtk(name)

    @test typeof(result) == Dict{Symbol, VecOrMat{<:Real}}
    @test result[:position] == position
    @test result[:Time] == [time]
    @test result[:Damage] == damage
    @test result[:Displacement] == displacement

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
        vtk["time", VTKFieldData()] = time
        vtk["damage", VTKPointData()] = damage
        vtk["random_point_data", VTKPointData()] = random_point_data
        vtk["random_field_data", VTKFieldData()] = random_field_data
    end
    result = read_vtk(name)

    @test typeof(result) == Dict{Symbol, VecOrMat{Float64}}
    @test result[:position] == position
    @test result[:time] == [time]
    @test result[:damage] == damage
    @test result[:random_point_data] == random_point_data
    @test result[:random_field_data] == random_field_data

    rm(name)
end

@testitem "wrong file type" begin
    @test_throws ArgumentError read_vtk("something.wrong")
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

@testitem "read pvtu file" begin
    using Peridynamics.WriteVTK

    root = joinpath(@__DIR__, "tempdir")
    isdir(root) && rm(root; recursive=true, force=true)
    bname = joinpath(root, "vtk_test")
    name = bname * ".pvtu"
    isfile(name) && rm(name)
    n_points = 10
    n_parts = 3
    position = [rand(3, n_points) for _ in 1:n_parts]
    cells = [Peridynamics.get_cells(n_points) for _ in 1:n_parts]
    damage = [rand(n_points) for _ in 1:n_parts]
    displacement = [rand(3, n_points) for _ in 1:n_parts]
    integer_test = [rand(1:100, n_points) for _ in 1:n_parts]
    time = rand() * 1000
    for i in 1:n_parts
        pvtk_grid(bname, position[i], cells[i]; part=i, nparts=n_parts) do vtk
            vtk["damage", VTKPointData()] = damage[i]
            vtk["displacement", VTKPointData()] = displacement[i]
            vtk["integer_test", VTKPointData()] = integer_test[i]
            vtk["time", VTKFieldData()] = time
        end
    end

    result = read_vtk(name)

    @test result[:position] ≈ reduce(hcat, position)
    @test result[:damage] ≈ reduce(vcat, damage)
    @test result[:displacement] ≈ reduce(hcat, displacement)
    @test result[:integer_test] ≈ reduce(vcat, integer_test)
    @test result[:time] ≈ [time]

    rm(root; recursive=true, force=true)
end
