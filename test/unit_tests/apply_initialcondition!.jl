using Peridynamics, Test

pc = PointCloud(1,1,1,0.25)
mat1 = BondBasedMaterial(; horizon = 1, rho = 1, E = 1, Gc = 1)
mat2 = ContinuumBasedMaterial(; horizon = 1, rho = 1, E = 1, nu = 1/4, Gc = 1)
body1 = Peridynamics.init_body(mat1, pc)
body2 = Peridynamics.init_body(mat2, pc)

vic1 = VelocityIC(2.0, 1:5, 1)

body1.velocity .= zeros(3, 64)
Peridynamics.apply_initialcondition!(body1, vic1)
@test body1.velocity[1, 1:5] == fill(2.0, 5)
