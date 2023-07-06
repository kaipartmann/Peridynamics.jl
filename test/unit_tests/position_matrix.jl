using Peridynamics, Test

let
    gridx = [1]
    gridy = [1]
    gridz = [1]
    position = Peridynamics.position_matrix(gridx, gridy, gridz)
    @test position == [1; 1; 1;;]
end

let
    gridx = [1, 2]
    gridy = [1, 2]
    gridz = [1, 2]
    position = Peridynamics.position_matrix(gridx, gridy, gridz)
    @test position == [
        1  2  1  2  1  2  1  2
        1  1  2  2  1  1  2  2
        1  1  1  1  2  2  2  2
    ]
end

let
    gridx = [1., 2.]
    gridy = [1., 2.]
    gridz = [1., 2.]
    position = Peridynamics.position_matrix(gridx, gridy, gridz)
    @test position == [
        1.0  2.0  1.0  2.0  1.0  2.0  1.0  2.0
        1.0  1.0  2.0  2.0  1.0  1.0  2.0  2.0
        1.0  1.0  1.0  1.0  2.0  2.0  2.0  2.0
    ]
end

let
    gridx = 1:2
    gridy = 1:2
    gridz = 1:2
    position = Peridynamics.position_matrix(gridx, gridy, gridz)
    @test position == [
        1  2  1  2  1  2  1  2
        1  1  2  2  1  1  2  2
        1  1  1  1  2  2  2  2
    ]
end

let
    gridx = 1.0:2.0
    gridy = 1.0:2.0
    gridz = 1.0:2.0
    position = Peridynamics.position_matrix(gridx, gridy, gridz)
    @test position == [
        1.0  2.0  1.0  2.0  1.0  2.0  1.0  2.0
        1.0  1.0  2.0  2.0  1.0  1.0  2.0  2.0
        1.0  1.0  1.0  1.0  2.0  2.0  2.0  2.0
    ]
end
