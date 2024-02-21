@testitem "find_points" begin
    #setup
    position = [
        0.0 1.0 0.0 0.0 1.0
        0.0 0.0 1.0 0.0 1.0
        0.0 0.0 0.0 1.0 1.0
    ]
    # find the points with different function names
    p = Peridynamics.find_points(p -> p[3] > 0.1 && p[2] > 0.1, position)
    @test p == [5]
    p = Peridynamics.find_points(p -> p[1] > 0.1 || p[2] > 0.1, position)
    @test p == [2, 3, 5]
    p = Peridynamics.find_points(p -> p[1] < -0.1 || p[2] < -0.1, position)
    @test p == Vector{Int}()
    px = Peridynamics.find_points(x -> x > 0.1, position)
    @test px == [2, 5]
    px = Peridynamics.find_points(x -> x < -0.1, position)
    @test px == Vector{Int}()
    py = Peridynamics.find_points(y -> y > 0.1, position)
    @test py == [3, 5]
    py = Peridynamics.find_points(y -> y < -0.1, position)
    @test py == Vector{Int}()
    pz = Peridynamics.find_points(z -> z > 0.1, position)
    @test pz == [4, 5]
    pz = Peridynamics.find_points(z -> z < -0.1, position)
    @test pz == Vector{Int}()

    @test_throws ArgumentError Peridynamics.find_points(k -> k > 0.1, position)
    @test_throws ArgumentError Peridynamics.find_points((k,x) -> x > 0.1, position)
    @test_throws ArgumentError Peridynamics.find_points(() -> 1, position)

    findfunc(x::Int) = x < 0
    findfunc(x::Float64) = x > 0
    @test_throws ArgumentError Peridynamics.find_points(findfunc, position)
end
