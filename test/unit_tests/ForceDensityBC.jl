using Peridynamics, Test

fdbc1 = ForceDensityBC(t -> 1.0, 1:5, 1)
@test fdbc1.fun(0) == 1.0
@test_throws MethodError fdbc1.fun() == 1.0
@test fdbc1.point_id_set == [1, 2, 3, 4, 5]
@test fdbc1.dim == 1
