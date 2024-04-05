@testitem "precrack!" begin
    # setup
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1, 1, 1, 1]
    mat = BBMaterial()
    body = Body(mat, position, volume)

    # test body creation
    @test body.n_points == 4
    @test body.position == position
    @test body.volume == volume
    @test body.fail_permit == fill(true, 4)

    # add point sets
    point_set!(body, :a, [1, 2])
    point_set!(body, :b, [2])
    point_set!(body, :c, [3])
    point_set!(body, :d, [4])

    # precrack! with dmg update = no bond filtering
    precrack!(body, :a, :c)
    @test body.point_sets_precracks == [Peridynamics.PointSetsPreCrack(:a, :c, false)]

    # precrack! without dmg update = bond filtering
    precrack!(body, :a, :d; update_dmg=false)
    @test body.point_sets_precracks == [Peridynamics.PointSetsPreCrack(:a, :c, false),
                                        Peridynamics.PointSetsPreCrack(:a, :d, true),]

    # precrack! with set :a and :b which intersect
    @test_throws ArgumentError precrack!(body, :a, :b)
end
