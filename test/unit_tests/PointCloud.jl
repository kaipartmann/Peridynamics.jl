using Peridynamics, Test

#-- PointCloud checks
#-- correct input
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_flag = BitVector([0, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1], "two" => [2])
    pc = PointCloud(n_points, position, volume, failure_flag, radius, point_sets)

    @test pc.n_points == n_points
    @test pc.position == position
    @test pc.volume == volume
    @test pc.failure_flag == failure_flag
    @test pc.radius == radius
    @test pc.point_sets == point_sets
end

#-- n_points smaller than zero
let
    n_points = 0
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_flag = BitVector([0, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Number of points `n_points` must be greater than zero!\n"
    @test_throws ErrorException(err_msg) PointCloud(n_points, position, volume,
                                                    failure_flag, radius, point_sets)
end

#-- position has wrong dimension
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
    ]
    volume = [1.0, 1.0]
    failure_flag = BitVector([0, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Incorrect dimensions of `position`!\n"
    err_msg *= "  Should be: (3, 2)\n"
    err_msg *= "  Evaluated: (2, 2)\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_flag, radius, point_sets)

    position = [
        1.0 2.0 2.1
        3.0 4.0 4.1
        5.0 6.0 6.1
    ]
    err_msg = "Incorrect dimensions of `position`!\n"
    err_msg *= "  Should be: (3, 2)\n"
    err_msg *= "  Evaluated: (3, 3)\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_flag, radius, point_sets)
end

#-- volume has wrong dimension
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0, 1.0]
    failure_flag = BitVector([0, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Incorrect length of `volume`!\n"
    err_msg *= "  Should be: 2\n"
    err_msg *= "  Evaluated: 3\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_flag, radius, point_sets)

    volume = [1.0]
    err_msg = "Incorrect length of `volume`!\n"
    err_msg *= "  Should be: 2\n"
    err_msg *= "  Evaluated: 1\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_flag, radius, point_sets)
end

#-- failure_flag has wrong dimension
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_flag = BitVector([0, 1, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Incorrect length of `failure_flag`!\n"
    err_msg *= "  Should be: 2\n"
    err_msg *= "  Evaluated: 3\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_flag, radius, point_sets)

    failure_flag = BitVector([0])
    err_msg = "Incorrect length of `failure_flag`!\n"
    err_msg *= "  Should be: 2\n"
    err_msg *= "  Evaluated: 1\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_flag, radius, point_sets)
end

#-- radius has wrong dimension
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_flag = BitVector([0, 1])
    radius = rand(3)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Incorrect length of `radius`!\n"
    err_msg *= "  Should be: 2\n"
    err_msg *= "  Evaluated: 3\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_flag, radius, point_sets)

    radius = rand(1)
    err_msg = "Incorrect length of `radius`!\n"
    err_msg *= "  Should be: 2\n"
    err_msg *= "  Evaluated: 1\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_flag, radius, point_sets)
end

#-- position contains NaN values
let
    n_points = 2
    position = [
        1.0 2.0
        NaN 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_flag = BitVector([0, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Matrix `position` contains NaN values!\n"
    @test_throws ErrorException(err_msg) PointCloud(n_points, position, volume,
                                                    failure_flag, radius, point_sets)
end

#-- volume contains NaN values
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, NaN]
    failure_flag = BitVector([0, 1])
    radius = fill((3 * 1.0 / (4π))^(1 / 3), 2)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Vector `volume` contains NaN values!\n"
    @test_throws ErrorException(err_msg) PointCloud(n_points, position, volume,
                                                    failure_flag, radius, point_sets)
end

#-- radius contains NaN values
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_flag = BitVector([0, 1])
    radius = [(3 * 1.0 / (4π))^(1 / 3), NaN]
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Vector `radius` contains NaN values!\n"
    @test_throws ErrorException(err_msg) PointCloud(n_points, position, volume,
                                                    failure_flag, radius, point_sets)
end

#-- BoundsError in poins_sets
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_flag = BitVector([0, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1, 3], "two" => [2])

    err_msg = "Invalid index [3] in point set `one`!\n"
    err_msg *= "Valid index range: 1:2\n"
    @test_throws ErrorException(err_msg) PointCloud(n_points, position, volume,
                                                    failure_flag, radius, point_sets)

    point_sets = Dict("one" => [1, 2], "two" => [-1, 2])

    err_msg = "Invalid index [-1] in point set `two`!\n"
    err_msg *= "Valid index range: 1:2\n"
    @test_throws ErrorException(err_msg) PointCloud(n_points, position, volume,
                                                    failure_flag, radius, point_sets)
end

#-- PointCloud from position and volume:
let
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    pc = PointCloud(position, volume)
    @test pc.n_points == 2
    @test pc.failure_flag == [true, true]
    @test pc.radius == Peridynamics.sphere_radius.(volume)
    @test isempty(pc.point_sets)
end

#-- PointCloud from length, width, height and Δx:
let
    pc = PointCloud(1, 1, 1, 0.5)
    @test pc.position == [
        -0.25  -0.25  -0.25  -0.25   0.25   0.25   0.25  0.25
        -0.25  -0.25   0.25   0.25  -0.25  -0.25   0.25  0.25
        -0.25   0.25  -0.25   0.25  -0.25   0.25  -0.25  0.25
    ]
    @test pc.n_points == 8
    @test pc.volume == fill(0.125, 8)
    @test pc.radius == Peridynamics.sphere_radius.(fill(0.125, 8))
    @test pc.failure_flag == fill(true, 8)
    @test isempty(pc.point_sets)
end

#-- PointCloud from length, width, height and Δx with center shift:
let
    pc = PointCloud(1, 1, 1, 0.5; center_x = 1, center_y = 1, center_z = 1)
    @test pc.position == [
        0.75  0.75  0.75  0.75  1.25  1.25  1.25  1.25
        0.75  0.75  1.25  1.25  0.75  0.75  1.25  1.25
        0.75  1.25  0.75  1.25  0.75  1.25  0.75  1.25
    ]
    @test pc.n_points == 8
    @test pc.volume == fill(0.125, 8)
    @test pc.radius == Peridynamics.sphere_radius.(fill(0.125, 8))
    @test pc.failure_flag == fill(true, 8)
    @test isempty(pc.point_sets)
end

#-- PointCloud from length, width, height and too big Δx:
let
    @test_throws ArgumentError PointCloud(1, 1, 0.1, 0.5)
end

#-- test the Base.show method of PointClouds
# for pc1 -> 2 points
let
    pc = PointCloud([1.0 2.0; 3.0 4.0; 5.0 6.0], [1.0, 1.0])
    io = IOBuffer()
    show(io, "text/plain", pc)
    msg = String(take!(io))
    @test msg == "2-points PointCloud" || msg == "2-points Peridynamics.PointCloud"
end

# for pc2 -> 8 points
let
    pc = PointCloud(1, 1, 1, 0.5)
    io = IOBuffer()
    show(io, "text/plain", pc)
    msg = String(take!(io))
    @test msg == "8-points PointCloud" || msg == "8-points Peridynamics.PointCloud"
end
