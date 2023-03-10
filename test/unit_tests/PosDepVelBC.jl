using Peridynamics, Test

pdvbc1 = PosDepVelBC((x, y, z, t) -> 1.0, 1:5, 1)
@test pdvbc1.fun(0, 0, 0, 0) == 1.0
@test_throws MethodError pdvbc1.fun() == 1.0
@test pdvbc1.point_id_set == [1, 2, 3, 4, 5]
@test pdvbc1.dim == 1
