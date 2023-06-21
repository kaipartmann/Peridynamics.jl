using Peridynamics, Test

fdbc1 = ForceDensityBC(t -> 1.0, 1:5, 1)
@test fdbc1.fun(0) == 1.0
@test_throws MethodError fdbc1.fun() == 1.0
@test fdbc1.point_id_set == [1, 2, 3, 4, 5]
@test fdbc1.dim == 1

io = IOBuffer()
show(io, "text/plain", fdbc1)
msg = String(take!(io))
@test msg == "5-points ForceDensityBC in x-direction (dim=1)" ||
      msg == "5-points Peridynamics.ForceDensityBC in x-direction (dim=1)"

fdbc2 = ForceDensityBC(t -> 1.0, 1:10, 2)
@test fdbc2.fun(0) == 1.0
@test_throws MethodError fdbc2.fun() == 1.0
@test fdbc2.point_id_set == [1:10;]
@test fdbc2.dim == 2

io = IOBuffer()
show(io, "text/plain", fdbc2)
msg = String(take!(io))
@test msg == "10-points ForceDensityBC in y-direction (dim=2)" ||
      msg == "10-points Peridynamics.ForceDensityBC in y-direction (dim=2)"

fdbc3 = ForceDensityBC(t -> 1.0, 1:20, 3)
@test fdbc3.fun(0) == 1.0
@test_throws MethodError fdbc3.fun() == 1.0
@test fdbc3.point_id_set == [1:20;]
@test fdbc3.dim == 3

io = IOBuffer()
show(io, "text/plain", fdbc3)
msg = String(take!(io))
@test msg == "20-points ForceDensityBC in z-direction (dim=3)" ||
      msg == "20-points Peridynamics.ForceDensityBC in z-direction (dim=3)" 