using Peridynamics, Test

pc = PointCloud(1,1,1,0.25)
mat1 = BBMaterial(; horizon = 1, rho = 1, E = 1, Gc = 1)
mat2 = CPDMaterial(; horizon = 1, rho = 1, E = 1, nu = 1/4, Gc = 1)
body1 = Peridynamics.init_body(mat1, pc)
body2 = Peridynamics.init_body(mat2, pc)

vbc1 = VelocityBC(t -> 1.0, 1:5, 1)
vbc2 = VelocityBC(t -> 2t, 1:5, 1)
pdvbc1 = PosDepVelBC((x, y, z, t) -> 1.0, 1:5, 1)
pdvbc2 = PosDepVelBC((x, y, z, t) -> x*y*z*t, 1:5, 2)
fdbc1 = ForceDensityBC(t -> 1.0, 1:5, 1)
fdbc2 = ForceDensityBC(t -> 2t, 1:5, 1)

t = 1.0

body1.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body1, vbc1, t)
@test body1.velocity_half[1, 1:5] == fill(1.0, 5)

body1.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body1, vbc2, t)
@test body1.velocity_half[1, 1:5] == fill(2.0, 5)

body1.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body1, pdvbc1, t)
@test body1.velocity_half[1, 1:5] == fill(1.0, 5)

body1.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body1, pdvbc2, t)
@test body1.velocity_half[2, 1:5] == [p[1] * p[2] * p[3] * t for p in eachcol(pc.position[:,1:5])]

body1.b_ext .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body1, fdbc1, t)
@test body1.b_ext[1, 1:5] == fill(1.0, 5)

body1.b_ext .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body1, fdbc2, t)
@test body1.b_ext[1, 1:5] == fill(2.0, 5)

body2.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body2, vbc1, t)
@test body2.velocity_half[1, 1:5] == fill(1.0, 5)

body2.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body2, vbc2, t)
@test body2.velocity_half[1, 1:5] == fill(2.0, 5)

body2.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body2, pdvbc1, t)
@test body2.velocity_half[1, 1:5] == fill(1.0, 5)

body2.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body2, pdvbc2, t)
@test body2.velocity_half[2, 1:5] == [p[1] * p[2] * p[3] * t for p in eachcol(pc.position[:,1:5])]

body2.b_ext .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body2, fdbc1, t)
@test body2.b_ext[1, 1:5] == fill(1.0, 5)

body2.b_ext .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body2, fdbc2, t)
@test body2.b_ext[1, 1:5] == fill(2.0, 5)

t = 2.0

body1.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body1, vbc1, t)
@test body1.velocity_half[1, 1:5] == fill(1.0, 5)

body1.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body1, vbc2, t)
@test body1.velocity_half[1, 1:5] == fill(4.0, 5)

body1.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body1, pdvbc1, t)
@test body1.velocity_half[1, 1:5] == fill(1.0, 5)

body1.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body1, pdvbc2, t)
@test body1.velocity_half[2, 1:5] == [p[1] * p[2] * p[3] * t for p in eachcol(pc.position[:,1:5])]

body1.b_ext .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body1, fdbc1, t)
@test body1.b_ext[1, 1:5] == fill(1.0, 5)

body1.b_ext .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body1, fdbc2, t)
@test body1.b_ext[1, 1:5] == fill(4.0, 5)

body2.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body2, vbc1, t)
@test body2.velocity_half[1, 1:5] == fill(1.0, 5)

body2.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body2, vbc2, t)
@test body2.velocity_half[1, 1:5] == fill(4.0, 5)

body2.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body2, pdvbc1, t)
@test body2.velocity_half[1, 1:5] == fill(1.0, 5)

body2.velocity_half .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body2, pdvbc2, t)
@test body2.velocity_half[2, 1:5] == [p[1] * p[2] * p[3] * t for p in eachcol(pc.position[:,1:5])]

body2.b_ext .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body2, fdbc1, t)
@test body2.b_ext[1, 1:5] == fill(1.0, 5)

body2.b_ext .= zeros(3, 64)
Peridynamics.apply_boundarycondition!(body2, fdbc2, t)
@test body2.b_ext[1, 1:5] == fill(4.0, 5)
