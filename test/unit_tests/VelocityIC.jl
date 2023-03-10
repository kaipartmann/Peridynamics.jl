using Peridynamics, Test

vic1 = VelocityIC(1.0, 1:5, 1)
@test vic1.val == 1.0
@test vic1.point_id_set == [1, 2, 3, 4, 5]
@test vic1.dim == 1
