using Peridynamics, Test
import Peridynamics: PointCloud

## PointCloud checks
# correct input
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_allowed = BitVector([0, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1], "two" => [2])
    pc = PointCloud(n_points, position, volume, failure_allowed, radius, point_sets)

    @test pc.n_points == n_points
    @test pc.position == position
    @test pc.volume == volume
    @test pc.failure_allowed == failure_allowed
    @test pc.radius == radius
    @test pc.point_sets == point_sets
end

## n_points smaller than zero
let
    n_points = 0
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_allowed = BitVector([0, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Number of points `n_points` must be greater than zero!\n"
    @test_throws ErrorException(err_msg) PointCloud(n_points, position, volume,
                                                    failure_allowed, radius, point_sets)
end

## position has wrong dimension
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
    ]
    volume = [1.0, 1.0]
    failure_allowed = BitVector([0, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Incorrect dimensions of `position`!\n"
    err_msg *= "  Should be: (3, 2)\n"
    err_msg *= "  Evaluated: (2, 2)\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_allowed, radius, point_sets)

    position = [
        1.0 2.0 2.1
        3.0 4.0 4.1
        5.0 6.0 6.1
    ]
    err_msg = "Incorrect dimensions of `position`!\n"
    err_msg *= "  Should be: (3, 2)\n"
    err_msg *= "  Evaluated: (3, 3)\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_allowed, radius, point_sets)
end

## volume has wrong dimension
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0, 1.0]
    failure_allowed = BitVector([0, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Incorrect length of `volume`!\n"
    err_msg *= "  Should be: 2\n"
    err_msg *= "  Evaluated: 3\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_allowed, radius, point_sets)

    volume = [1.0]
    err_msg = "Incorrect length of `volume`!\n"
    err_msg *= "  Should be: 2\n"
    err_msg *= "  Evaluated: 1\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_allowed, radius, point_sets)
end

## failure_allowed has wrong dimension
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_allowed = BitVector([0, 1, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Incorrect length of `failure_allowed`!\n"
    err_msg *= "  Should be: 2\n"
    err_msg *= "  Evaluated: 3\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_allowed, radius, point_sets)

    failure_allowed = BitVector([0])
    err_msg = "Incorrect length of `failure_allowed`!\n"
    err_msg *= "  Should be: 2\n"
    err_msg *= "  Evaluated: 1\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_allowed, radius, point_sets)
end

## radius has wrong dimension
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_allowed = BitVector([0, 1])
    radius = rand(3)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Incorrect length of `radius`!\n"
    err_msg *= "  Should be: 2\n"
    err_msg *= "  Evaluated: 3\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_allowed, radius, point_sets)

    radius = rand(1)
    err_msg = "Incorrect length of `radius`!\n"
    err_msg *= "  Should be: 2\n"
    err_msg *= "  Evaluated: 1\n"
    @test_throws DimensionMismatch(err_msg) PointCloud(n_points, position, volume,
                                                       failure_allowed, radius, point_sets)
end

## position contains NaN values
let
    n_points = 2
    position = [
        1.0 2.0
        NaN 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_allowed = BitVector([0, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Matrix `position` contains NaN values!\n"
    @test_throws ErrorException(err_msg) PointCloud(n_points, position, volume,
                                                    failure_allowed, radius, point_sets)
end

## volume contains NaN values
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, NaN]
    failure_allowed = BitVector([0, 1])
    radius = fill((3 * 1.0 / (4π))^(1 / 3), 2)
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Vector `volume` contains NaN values!\n"
    @test_throws ErrorException(err_msg) PointCloud(n_points, position, volume,
                                                    failure_allowed, radius, point_sets)
end

## radius contains NaN values
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_allowed = BitVector([0, 1])
    radius = [(3 * 1.0 / (4π))^(1 / 3), NaN]
    point_sets = Dict("one" => [1], "two" => [2])

    err_msg = "Vector `radius` contains NaN values!\n"
    @test_throws ErrorException(err_msg) PointCloud(n_points, position, volume,
                                                    failure_allowed, radius, point_sets)
end

## BoundsError in poins_sets
let
    n_points = 2
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    failure_allowed = BitVector([0, 1])
    radius = Peridynamics.sphere_radius.(volume)
    point_sets = Dict("one" => [1, 3], "two" => [2])

    err_msg = "Invalid index [3] in point set `one`!\n"
    err_msg *= "Valid index range: 1:2\n"
    @test_throws ErrorException(err_msg) PointCloud(n_points, position, volume,
                                                    failure_allowed, radius, point_sets)

    point_sets = Dict("one" => [1, 2], "two" => [-1, 2])

    err_msg = "Invalid index [-1] in point set `two`!\n"
    err_msg *= "Valid index range: 1:2\n"
    @test_throws ErrorException(err_msg) PointCloud(n_points, position, volume,
                                                    failure_allowed, radius, point_sets)
end

## PointCloud from position and volume:
let
    position = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume = [1.0, 1.0]
    pc = PointCloud(position, volume)
    @test pc.n_points == 2
    @test pc.failure_allowed == [true, true]
    @test pc.radius == Peridynamics.sphere_radius.(volume)
    @test isempty(pc.point_sets)
end

## PointCloud from length, width, height and Δx:
let
    pc = PointCloud(1, 1, 1, 0.5)
    @test pc.position == [
        -0.25   0.25  -0.25   0.25  -0.25   0.25  -0.25  0.25
        -0.25  -0.25   0.25   0.25  -0.25  -0.25   0.25  0.25
        -0.25  -0.25  -0.25  -0.25   0.25   0.25   0.25  0.25
    ]
    @test pc.n_points == 8
    @test pc.volume == fill(0.125, 8)
    @test pc.radius == Peridynamics.sphere_radius.(fill(0.125, 8))
    @test pc.failure_allowed == fill(true, 8)
    @test isempty(pc.point_sets)
end

## PointCloud from length, width, height and Δx with center shift:
let
    pc = PointCloud(1, 1, 1, 0.5; center_x = 1, center_y = 1, center_z = 1)
    @test pc.position == [
        0.75  1.25  0.75  1.25  0.75  1.25  0.75  1.25
        0.75  0.75  1.25  1.25  0.75  0.75  1.25  1.25
        0.75  0.75  0.75  0.75  1.25  1.25  1.25  1.25
    ]
    @test pc.n_points == 8
    @test pc.volume == fill(0.125, 8)
    @test pc.radius == Peridynamics.sphere_radius.(fill(0.125, 8))
    @test pc.failure_allowed == fill(true, 8)
    @test isempty(pc.point_sets)
end

## PointCloud from length, width, height and too big Δx:
let
    @test_throws ArgumentError PointCloud(1, 1, 0.1, 0.5)
end

## test the Base.show method of PointClouds
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

## pcmerge
let
    position1 = [
        -0.25  -0.25  -0.25  -0.25   0.25   0.25   0.25  0.25
        -0.25  -0.25   0.25   0.25  -0.25  -0.25   0.25  0.25
        -0.25   0.25  -0.25   0.25  -0.25   0.25  -0.25  0.25
    ]
    volume1 = fill(1.0, 8)
    point_sets1 = Dict("set11" => [1,2], "set12" => [3,4])
    pc1 = PointCloud(position1, volume1, point_sets1)
    pc1.failure_allowed[1:2] .= false
    @test pc1.position == position1
    @test pc1.n_points == 8
    @test pc1.volume == volume1
    @test pc1.radius == Peridynamics.sphere_radius.(volume1)
    @test pc1.failure_allowed == BitVector([0,0,1,1,1,1,1,1])
    @test pc1.point_sets == point_sets1

    position1_set11 = position1[:, point_sets1["set11"]]
    position1_set12 = position1[:, point_sets1["set12"]]
    @test pc1.position[:, pc1.point_sets["set11"]] ≈ position1_set11
    @test pc1.position[:, pc1.point_sets["set12"]] ≈ position1_set12

    position2 = [
        0.75   0.75   0.75   0.75   1.25   1.25   1.25  1.25
        -0.25  -0.25   0.25   0.25  -0.25  -0.25   0.25  0.25
        -0.25   0.25  -0.25   0.25  -0.25   0.25  -0.25  0.25
    ]
    volume2 = fill(1.0, 8)
    point_sets2 = Dict("set21" => [5,6], "set22" => [7,8])
    pc2 = PointCloud(position2, volume2, point_sets2)
    pc2.failure_allowed[7:8] .= false
    @test pc2.position == position2
    @test pc2.n_points == 8
    @test pc2.volume == volume2
    @test pc2.radius == Peridynamics.sphere_radius.(volume2)
    @test pc2.failure_allowed == BitVector([1,1,1,1,1,1,0,0])
    @test pc2.point_sets == point_sets2

    position2_set21 = position2[:, point_sets2["set21"]]
    position2_set22 = position2[:, point_sets2["set22"]]
    @test pc2.position[:, pc2.point_sets["set21"]] ≈ position2_set21
    @test pc2.position[:, pc2.point_sets["set22"]] ≈ position2_set22

    pc = Peridynamics.pcmerge([pc1, pc2])
    @test pc.n_points == 16
    @test pc.position == hcat(position1, position2)
    @test pc.volume == vcat(volume1, volume2)
    @test pc.failure_allowed == [0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0]
    @test pc.radius == Peridynamics.sphere_radius.(vcat(volume1, volume2))
    @test pc.point_sets == Dict("set11" => [1,2], "set12" => [3,4], "set21" => [13,14],
                                "set22" => [15,16])
end
