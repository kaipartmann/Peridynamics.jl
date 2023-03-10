using Peridynamics, Test

position1 = [
    -0.25  -0.25  -0.25  -0.25   0.25   0.25   0.25  0.25
    -0.25  -0.25   0.25   0.25  -0.25  -0.25   0.25  0.25
    -0.25   0.25  -0.25   0.25  -0.25   0.25  -0.25  0.25
]
volume1 = fill(1.0, 8)
point_sets1 = Dict("set11" => [1,2], "set12" => [3,4])
pc1 = PointCloud(position1, volume1, point_sets1)
pc1.failure_flag[1:2] .= false
@test pc1.position == position1
@test pc1.n_points == 8
@test pc1.volume == volume1
@test pc1.radius == Peridynamics.sphere_radius.(volume1)
@test pc1.failure_flag == BitVector([0,0,1,1,1,1,1,1])
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
pc2.failure_flag[7:8] .= false
@test pc2.position == position2
@test pc2.n_points == 8
@test pc2.volume == volume2
@test pc2.radius == Peridynamics.sphere_radius.(volume2)
@test pc2.failure_flag == BitVector([1,1,1,1,1,1,0,0])
@test pc2.point_sets == point_sets2

position2_set21 = position2[:, point_sets2["set21"]]
position2_set22 = position2[:, point_sets2["set22"]]
@test pc2.position[:, pc2.point_sets["set21"]] ≈ position2_set21
@test pc2.position[:, pc2.point_sets["set22"]] ≈ position2_set22

pc = pcmerge([pc1, pc2])
@test pc.n_points == 16
@test pc.position == hcat(position1, position2)
@test pc.volume == vcat(volume1, volume2)
@test pc.failure_flag == [0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0]
@test pc.radius == Peridynamics.sphere_radius.(vcat(volume1, volume2))
@test pc.point_sets == Dict("set11" => [1,2], "set12" => [3,4], "set21" => [13,14],
                            "set22" => [15,16])
