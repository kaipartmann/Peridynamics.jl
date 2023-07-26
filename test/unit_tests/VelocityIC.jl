using Peridynamics, Test

vic1 = VelocityIC(1.0, 1:5, 1)
@test vic1.val == 1.0
@test vic1.point_id_set == [1, 2, 3, 4, 5]
@test vic1.dim == 1

io = IOBuffer()
show(io, "text/plain", vic1)
msg = String(take!(io))
@test msg == "5-points VelocityIC in x-direction (dim=1)" || 
      msg == "5-points Peridynamics.VelocityIC in x-direction (dim=1)"

vic2 = VelocityIC(1.0, 1:10, 2)
@test vic2.val == 1.0
@test vic2.point_id_set == [1:10;]
@test vic2.dim == 2

io = IOBuffer()
show(io, "text/plain", vic2)
msg = String(take!(io))
@test msg == "10-points VelocityIC in y-direction (dim=2)" ||
      msg == "10-points Peridynamics.VelocityIC in y-direction (dim=2)"

vic3 = VelocityIC(1.0, 1:20, 3)
@test vic3.val == 1.0
@test vic3.point_id_set == [1:20;]
@test vic3.dim == 3

io = IOBuffer()
show(io, "text/plain", vic3)
msg = String(take!(io))
@test msg == "20-points VelocityIC in z-direction (dim=3)" || 
      msg == "20-points Peridynamics.VelocityIC in z-direction (dim=3)"