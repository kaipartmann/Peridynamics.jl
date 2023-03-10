using Peridynamics, Test

vbc1 = VelocityBC(t -> 1.0, 1:5, 1)
@test vbc1.fun(0) == 1.0
@test_throws MethodError vbc1.fun() == 1.0
@test vbc1.point_id_set == [1, 2, 3, 4, 5]
@test vbc1.dim == 1
