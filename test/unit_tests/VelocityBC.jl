using Peridynamics, Test

vbc1 = VelocityBC(t -> 1.0, 1:5, 1)
@test vbc1.fun(0) == 1.0
@test_throws MethodError vbc1.fun() == 1.0
@test vbc1.point_id_set == [1, 2, 3, 4, 5]
@test vbc1.dim == 1

io = IOBuffer()
show(io, "text/plain", vbc1)
msg = String(take!(io))
@test msg == "5-points VelocityBC in x-direction (dim=1)" || 
      msg == "5-points Peridynamics.VelocityBC in x-direction (dim=1)"