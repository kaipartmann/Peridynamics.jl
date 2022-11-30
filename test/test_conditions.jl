using Peridynamics
using Test

#------------------------------------------------------------------------------------------
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
