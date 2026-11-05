@testitem "get_points" begin
    import Peridynamics.AbaqusMeshConverter: get_points

    nodes1 = Dict(
        1 => [0.0, 0.0, 0.0],
        2 => [1.0, 0.0, 0.0],
        3 => [0.0, 1.0, 0.0],
        4 => [0.0, 0.0, 1.0],
    )
    elements1 = Dict(
        1 => [1, 2, 3, 4]
    )
    element_types1 = Dict(
        1 => :Tet4
    )
    midpoints1, volumes1 = get_points(nodes1, elements1, element_types1)
    @test midpoints1 == [0.25; 0.25; 0.25;;]
    @test volumes1 == [1/6]

    wrong_elem1 = Dict(
        1 => [1, 2, 3, 4, 5]
    )
    @test_throws DimensionMismatch get_points(nodes1, wrong_elem1, element_types1)

    nodes2 = Dict(
        1 => [0.0, 1.0, 1.0],
        2 => [0.0, 0.0, 1.0],
        3 => [0.0, 0.0, 0.0],
        4 => [0.0, 1.0, 0.0],
        5 => [1.0, 1.0, 1.0],
        6 => [1.0, 0.0, 1.0],
        7 => [1.0, 0.0, 0.0],
        8 => [1.0, 1.0, 0.0],
    )
    elements2 = Dict(
        1 => [1, 2, 3, 4, 5, 6, 7, 8]
    )
    element_types2 = Dict(
        1 => :Hex8
    )
    midpoints2, volumes2 = get_points(nodes2, elements2, element_types2)
    @test midpoints2 == [0.5; 0.5; 0.5;;]
    @test volumes2 == [1]

    wrong_elem2 = Dict(
        1 => [1, 2, 3, 4, 5, 6, 7]
    )
    @test_throws DimensionMismatch get_points(nodes1, wrong_elem1, element_types1)

    wrong_elemtype = Dict(
        1 => :Tet8
    )
    @test_throws DomainError get_points(nodes1, elements1, wrong_elemtype)
end

@testitem "midpoint" begin
    import Peridynamics.AbaqusMeshConverter: midpoint

    @test midpoint([1], [2], [3], [4]) == [2.5]
    @test midpoint(([0, 0, 0] for _ in 1:4)...) == [0, 0, 0]
    @test midpoint(([0, 0, 0] for _ in 1:8)...) == [0, 0, 0]
    @test midpoint(([i, i+1, i+2] for i in 1:4)...) == [2.5, 3.5, 4.5]
    @test midpoint(([i, i+1, i+2] for i in 1:8)...) == [4.5, 5.5, 6.5]
    @test_throws MethodError midpoint([1], [2], [3])
end

@testitem "tetvol" begin
    import Peridynamics.AbaqusMeshConverter: tetvol

    @test tetvol([0,0,0], [1,0,0], [0,1,0], [0,0,1]) == 1/6
    @test tetvol(([0,0,0] for _ in 1:4)...) == 0
end

@testitem "read_inp CubeC3D8" begin
    file = joinpath(@__DIR__, "models", "CubeC3D8.inp")
    pos, vol, sets = read_inp(file)
    @test size(pos) == (3, 125)
    @test length(vol) == 125
    @test vol == fill(4^3, 125)
    @test maximum(pos[1,:]) == 18
    @test maximum(pos[2,:]) == 18
    @test maximum(pos[3,:]) == 18
    @test minimum(pos[3,:]) == 2
    @test minimum(pos[3,:]) == 2
    @test minimum(pos[3,:]) == 2
    @test sets["l"] == 101:125
    @test sets["r"] == 1:25
end

@testitem "read_inp CubeC3D4" begin
    file = joinpath(@__DIR__, "models", "CubeC3D4.inp")
    pos, vol, sets = read_inp(file)
    @test size(pos) == (3, 1050)
    @test length(vol) == 1050
    @test sum((vol .> 0)) == 1050
    @test maximum(vol) ≈ 22.55287594750887
    @test minimum(vol) ≈ 2.9106088794695766
    @test maximum(pos) ≈ 19.55310155
    @test minimum(pos) ≈ 0.416666805
end

@testitem "read_inp wrong input" begin
    @test_throws ArgumentError read_inp("something.wrong")

    file_with_unknown_type = joinpath(@__DIR__, "models", "CubeC3D10.inp")
    @test_throws DomainError read_inp(file_with_unknown_type)
end
