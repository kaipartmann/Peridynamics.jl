using Peridynamics, Test

pc1 = read_inp(joinpath(@__DIR__, "models", "CubeC3D8.inp"))
@test pc1.n_points == 125
@test pc1.volume == fill(4^3, 125)
@test maximum(pc1.position[1,:]) == 18
@test maximum(pc1.position[2,:]) == 18
@test maximum(pc1.position[3,:]) == 18
@test minimum(pc1.position[3,:]) == 2
@test minimum(pc1.position[3,:]) == 2
@test minimum(pc1.position[3,:]) == 2
@test pc1.point_sets["l"] == 101:125
@test pc1.point_sets["r"] == 1:25

pc2 = read_inp(joinpath(@__DIR__, "models", "CubeC3D4.inp"))
@test pc2.n_points == 1050
@test sum((pc2.volume .> 0)) == 1050
@test maximum(pc2.volume) ≈ 22.55287594750887
@test minimum(pc2.volume) ≈ 2.9106088794695766
@test maximum(pc2.position) ≈ 19.55310155
@test minimum(pc2.position) ≈ 0.416666805

@test_throws AssertionError read_inp("something.wrong")
@test_throws DomainError read_inp(joinpath(@__DIR__, "models", "CubeC3D10.inp"))
