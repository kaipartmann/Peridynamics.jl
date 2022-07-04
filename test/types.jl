using Peridynamics
using Test

##------------------------------------------------------------------------------------------
# PointCloud

position1 = [
    1.0 2.0
    3.0 4.0
    5.0 6.0
]
volume1 = [1.0, 1.0]
pc1 = PointCloud(position1, volume1)
@test pc1.n_points == 2
@test pc1.failure_flag == [true, true]
@test pc1.radius == Peridynamics.sphere_radius.(volume1)
@test isempty(pc1.point_sets)
@test_throws DimensionMismatch PointCloud([1.0 2.0; 3.0 4.0], volume1)

pc2 = PointCloud(1, 1, 1, 0.5)
@test pc2.position == [
    -0.25  -0.25  -0.25  -0.25   0.25   0.25   0.25  0.25
    -0.25  -0.25   0.25   0.25  -0.25  -0.25   0.25  0.25
    -0.25   0.25  -0.25   0.25  -0.25   0.25  -0.25  0.25
]
@test pc2.n_points == 8
@test pc2.volume == fill(0.125, 8)
@test pc2.radius == Peridynamics.sphere_radius.(fill(0.125, 8))
@test pc2.failure_flag == fill(true, 8)
@test isempty(pc2.point_sets)

pc3 = PointCloud(1, 1, 1, 0.5; center_x = 1, center_y = 1, center_z = 1)
@test pc3.position == [
    0.75  0.75  0.75  0.75  1.25  1.25  1.25  1.25
    0.75  0.75  1.25  1.25  0.75  0.75  1.25  1.25
    0.75  1.25  0.75  1.25  0.75  1.25  0.75  1.25
]
@test pc3.n_points == 8
@test pc3.volume == fill(0.125, 8)
@test pc3.radius == Peridynamics.sphere_radius.(fill(0.125, 8))
@test pc3.failure_flag == fill(true, 8)
@test isempty(pc3.point_sets)

##------------------------------------------------------------------------------------------
# TimeDiscretization

td1 = TimeDiscretization(1)
@test td1.n_timesteps == 1
@test td1.Δt == -1
@test td1.alg == :verlet

td2 = TimeDiscretization(2; alg=:dynrelax)
@test td2.n_timesteps == 2
@test td2.Δt == 1
@test td2.alg == :dynrelax

@test_throws DomainError TimeDiscretization(2; alg = :somethingwrong)

td3 = TimeDiscretization(3, 0.1)
@test td3.n_timesteps == 3
@test td3.Δt == 0.1
@test td3.alg == :verlet

msg = "alg = :dynrelax -> will skip given value of Δt! Check your input values!"
td4 = @test_warn msg TimeDiscretization(3, 0.2; alg=:dynrelax)
@test td4.Δt == 1
@test td4.alg == :dynrelax

@test_throws DomainError TimeDiscretization(4, 0.4; alg = :somethingwrong)


##------------------------------------------------------------------------------------------
# ExportSettings

es1 = ExportSettings()
@test es1.path == ""
@test es1.exportfreq == 0
@test es1.resultfile_prefix == ""
@test es1.logfile == ""
@test es1.exportflag == false

es2 = ExportSettings("test/path", 10)
@test es2.path == "test/path"
@test es2.exportfreq == 10
@test es2.resultfile_prefix == ""
@test es2.logfile == ""
@test es2.exportflag == true


##------------------------------------------------------------------------------------------
# Boundary Conditions and Initial Conditions

vbc1 = VelocityBC(t -> 1.0, 1:5, 1)
@test vbc1.fun(0) == 1.0
@test_throws MethodError vbc1.fun() == 1.0
@test vbc1.point_id_set == [1, 2, 3, 4, 5]
@test vbc1.dim == 1

vic1 = VelocityIC(1.0, 1:5, 1)
@test vic1.val == 1.0
@test vic1.point_id_set == [1, 2, 3, 4, 5]
@test vic1.dim == 1

fdbc1 = ForceDensityBC(t -> 1.0, 1:5, 1)
@test fdbc1.fun(0) == 1.0
@test_throws MethodError fdbc1.fun() == 1.0
@test fdbc1.point_id_set == [1, 2, 3, 4, 5]
@test fdbc1.dim == 1

pdvbc1 = PosDepVelBC((x, y, z, t) -> 1.0, 1:5, 1)
@test pdvbc1.fun(0, 0, 0, 0) == 1.0
@test_throws MethodError pdvbc1.fun() == 1.0
@test pdvbc1.point_id_set == [1, 2, 3, 4, 5]
@test pdvbc1.dim == 1
