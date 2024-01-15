using Peridynamics, Test

mat1 = CPDMaterial(; horizon = 1, rho = 1, E = 1, nu = 1/4, Gc = 1)
@test mat1.δ == 1
@test mat1.rho == 1
@test mat1.E == 1
@test mat1.nu == 1/4
@test mat1.G ≈ 1 / (2 * (1 + 1/4))
@test mat1.K ≈ 1 / (3 * (1 - 2 * 1/4))
@test mat1.C1 == 30/π * 1 / (2 * (1 + 1/4))/1^4
@test mat1.C2 == 0
@test mat1.C3 == 0
@test mat1.Gc == 1
@test mat1.εc == sqrt(5.0 * 1 / (9.0 * 1 / (3 * (1 - 2 * 1/4)) * 1))

mat2 = CPDMaterial(; horizon = 0, rho = 0, E = 0, nu = 0, Gc = 0)
@test mat2.δ == 0
@test mat2.rho == 0
@test mat2.E == 0
@test mat2.nu == 0
@test mat2.G == 0
@test mat2.K == 0
@test isnan(mat2.C1)
@test mat2.C2 == 0
@test isnan(mat2.C3)
@test mat2.Gc == 0
@test isnan(mat2.εc)

mat3 = CPDMaterial(; horizon = 1, rho = 8e-6, E = 210e3, nu = 0.3, epsilon_c = 0.01)
@test mat3.δ == 1
@test mat3.rho == 8e-6
@test mat3.E == 210e3
@test mat3.nu == 0.3
@test mat3.G ≈ 210e3 / (2 * (1 + 0.3))
@test mat3.K ≈ 210e3 / (3 * (1 - 2 * 0.3))
@test mat3.C1 == 30/π * 210e3 / (2 * (1 + 0.3))/1^4
@test mat3.C2 == 0
@test mat3.C3 == 32/π^4 * (210e3 * 0.3 /((1 + 0.3) * (1 - 2*0.3)) - 210e3 / (2 * (1 + 0.3)))/1^12
@test mat3.Gc ≈ 9.0/5.0 * 210e3 / (3 * (1 - 2 * 0.3)) * 1 * 0.01^2
@test mat3.εc == 0.01

err = ArgumentError("Duplicate definition: define either Gc or epsilon_c, not both!")
@test_throws err CPDMaterial(; horizon = 1, rho = 8e-6, E = 210e3, nu=0.3, epsilon_c = 0.01, Gc = 1)

err = ArgumentError("Either Gc or epsilon_c have to be defined!")
@test_throws err CPDMaterial(; horizon = 1, rho = 8e-6, E = 210e3, nu=0.3)

mat4 = CPDMaterial(horizon=1.5, rho=7850.0, E=210e9, nu=0.3, Gc=1.0)
io = IOBuffer()
show(io, "text/plain", mat4)
msg = String(take!(io))

correct_msg1 = "CPDMaterial:\n  δ: 1.5\n  rho: 7850\n  E: 2.1e+11\n  nu: 0.3\n  G: 8.07692e+10\n  K: 1.75e+11\n  C1: 1.52353e+11\n  C2: 0\n  C3: 1.02252e+08\n  Gc: 1\n  εc: 1.45479e-06\n"
correct_msg2 = "Peridynamics.CPDMaterial:\n  δ: 1.5\n  rho: 7850\n  E: 2.1e+11\n  nu: 0.3\n  G: 8.07692e+10\n  K: 1.75e+11\n  C1: 1.52353e+11\n  C2: 0\n  C3: 1.02252e+08\n  Gc: 1\n  εc: 1.45479e-06\n"

@test msg == correct_msg1 || msg == correct_msg2

warn_msg = "CPD parameters choosen manually!\nBe careful when adjusting CPD parameters to avoid unexpected outcomes!"
mat5 = @test_logs (:warn, warn_msg) CPDMaterial(; horizon = 1, rho = 1, E = 1, nu = 1/4, Gc = 1, C1=1)
@test mat5.δ == 1
@test mat5.rho == 1
@test mat5.E == 1
@test mat5.nu == 1/4
@test mat5.G ≈ 1 / (2 * (1 + 1/4))
@test mat5.K ≈ 1 / (3 * (1 - 2 * 1/4))
@test mat5.C1 == 1
@test mat5.C2 == 0
@test mat5.C3 == 0
@test mat5.Gc == 1
@test mat5.εc == sqrt(5.0 * 1 / (9.0 * 1 / (3 * (1 - 2 * 1/4)) * 1))
