using Peridynamics, Test

pc = PointCloud(1, 0.5, 0.5, 0.25)
pc.failure_flag .= false
mat1 = CPDMaterial(; horizon=0.8, rho=8000, E=210e9, nu=0.2, Gc=1, C1=1, C2=1, C3=1)
mat2 = deepcopy(mat1)
mat_same = MultiMaterial((mat1, mat2), ifelse.(pc.position[1,:] .< 0, 1, 2))

body1 = Peridynamics.init_body(mat1, pc)
body2 = Peridynamics.init_body(mat_same, pc)

@test body1.n_bonds == body2.n_bonds
@test body1.n_two_ni == body2.n_two_ni
@test body1.n_three_ni == body2.n_three_ni

mat3 = CPDMaterial(; horizon=0.8, rho=8000, E=210e9, nu=0.2, Gc=1, C1=1, C2=0, C3=0)
mat_half = MultiMaterial((mat1, mat3), ifelse.(pc.position[1,:] .< 0, 1, 2))
body3 = Peridynamics.init_body(mat1, pc)
body4 = Peridynamics.init_body(mat_half, pc)

@test body3.n_bonds == body4.n_bonds
@test body3.n_two_ni == 2 * body4.n_two_ni
@test body3.n_three_ni == 2 * body4.n_three_ni

fun1(t) = 0.0001
fun2(t) = -0.0001
bc1 = VelocityBC(fun1,1:4,1)
bc2 = VelocityBC(fun2,13:16,1)
bcs = [bc1, bc2]

td = DynamicRelaxation(100)
es = ExportSettings()

job = deepcopy(PDSingleBodyAnalysis(name="mmi_1",pc=pc,mat=mat1,bcs=bcs,td=td,es=es))
job_mm = deepcopy(PDSingleBodyAnalysis(name="mmi_same",pc=pc,mat=mat_same,bcs=bcs,td=td,es=es))

body = submit(job)
body_mm = submit(job_mm)

@test body.n_bonds == body_mm.n_bonds
@test body.n_two_ni == body_mm.n_two_ni
@test body.n_three_ni == body_mm.n_three_ni
