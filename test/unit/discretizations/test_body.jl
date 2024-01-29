using Peridynamics, Test

##
let
    # setup
    n_points = 10
    position, volume = rand(3, n_points), rand(n_points)
    body = Body(position, volume)

    # test body creation
    @test body.n_points == n_points
    @test body.position == position
    @test body.volume == volume
    @test body.failure_allowed == BitVector(fill(true, n_points))
    @test isa(body.psh, Peridynamics.UndefPointSetHandler)
end

##
let
    # setup
    position = [
        0.0 1.0 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 0.0 0.0 1.0
    ]
    volume = [1, 1, 1, 1]
    body = Body(position, volume)

    # test body creation
    @test body.n_points == 4
    @test body.position == position
    @test body.volume == volume
    @test body.failure_allowed == BitVector(fill(true, 4))
    @test isa(body.psh, Peridynamics.UndefPointSetHandler)

    # add point set
    point_set!(body, :a, 1:2)
    @test isa(body.psh, Peridynamics.NoMatPointSetHandler)
    @test body.psh.point_sets == Dict(:a => 1:2)

    # add another point set via function definition
    point_set!(x -> x > 0.5, body, :b)
    @test isa(body.psh, Peridynamics.NoMatPointSetHandler)
    @test body.psh.point_sets == Dict(:a => 1:2, :b => [2])

    # add point set with do syntax
    point_set!(body, :c) do p
        p[3] > 0.0
    end
    @test isa(body.psh, Peridynamics.NoMatPointSetHandler)
    @test body.psh.point_sets == Dict(:a => 1:2, :b => [2], :c => [4])

    # point_set!
    @test_throws BoundsError point_set!(body, :d, 1:5)
end
