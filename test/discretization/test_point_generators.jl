@testitem "uniform_box" begin
    pos, vol = uniform_box(1, 1, 1, 0.5)
    @test pos == [
        -0.25   0.25  -0.25   0.25  -0.25   0.25  -0.25  0.25
        -0.25  -0.25   0.25   0.25  -0.25  -0.25   0.25  0.25
        -0.25  -0.25  -0.25  -0.25   0.25   0.25   0.25  0.25
    ]
    @test vol == fill(0.125, 8)

    pos, vol = uniform_box(1, 1, 1, 0.5; center_x=1, center_y=1, center_z=1)
    @test pos == [
        0.75  1.25  0.75  1.25  0.75  1.25  0.75  1.25
        0.75  0.75  1.25  1.25  0.75  0.75  1.25  1.25
        0.75  0.75  0.75  0.75  1.25  1.25  1.25  1.25
    ]
    @test vol == fill(0.125, 8)
end

@testitem "uniform_sphere" begin
    pos, vol = uniform_sphere(1, 0.5)
    @test pos == [
        -0.25   0.25  -0.25   0.25  -0.25   0.25  -0.25  0.25
        -0.25  -0.25   0.25   0.25  -0.25  -0.25   0.25  0.25
        -0.25  -0.25  -0.25  -0.25   0.25   0.25   0.25  0.25
    ]
    @test vol == fill(0.125, 8)

    pos, vol = uniform_sphere(1, 0.5; center_x=1, center_y=1, center_z=1)
    @test pos == [
        0.75  1.25  0.75  1.25  0.75  1.25  0.75  1.25
        0.75  0.75  1.25  1.25  0.75  0.75  1.25  1.25
        0.75  0.75  0.75  0.75  1.25  1.25  1.25  1.25
    ]
    @test vol == fill(0.125, 8)
end
