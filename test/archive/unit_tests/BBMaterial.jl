using Peridynamics, Test

mat1 = BBMaterial(; horizon = 1, rho = 1, E = 1, Gc = 1)
@test mat1.δ == 1
@test mat1.rho == 1
@test mat1.E == 1
@test mat1.nu == 1/4
@test mat1.G == 1 / (2 * (1 + 1/4))
@test mat1.K ≈ 1 / (3 * (1 - 2 * 1/4))
@test mat1.bc ≈ 18 * 1 / (3 * (1 - 2 * 1/4)) / (π * 1^4)
@test mat1.Gc == 1
@test mat1.εc == sqrt(5.0 * 1 / (9.0 * 1 / (3 * (1 - 2 * 1/4)) * 1))

mat2 = BBMaterial(; horizon = 0, rho = 0, E = 0, Gc = 0)
@test mat2.δ == 0
@test mat2.rho == 0
@test mat2.E == 0
@test mat2.nu == 1/4
@test mat2.G == 0
@test mat2.K == 0
@test isnan(mat2.bc)
@test mat2.Gc == 0
@test isnan(mat2.εc)

mat3 = BBMaterial(; horizon = 1, rho = 8e-6, E = 210e3, epsilon_c = 0.01)
@test mat3.δ == 1
@test mat3.rho == 8e-6
@test mat3.E == 210e3
@test mat3.nu == 1/4
@test mat3.G ≈ 210e3 / (2 * (1 + 0.25))
@test mat3.K ≈ 210e3 / (3 * (1 - 2 * 0.25))
@test mat3.bc ≈ 18 * 210e3 / (3 * (1 - 2 * 0.25)) / (π * 1^4)
@test mat3.Gc ≈ 9.0/5.0 * 210e3 / (3 * (1 - 2 * 0.25)) * 1 * 0.01^2
@test mat3.εc == 0.01

err = ArgumentError("Duplicate definition: define either Gc or epsilon_c, not both!")
@test_throws err BBMaterial(; horizon = 1, rho = 8e-6, E = 210e3, epsilon_c = 0.01, Gc = 1)

err = ArgumentError("Either Gc or epsilon_c have to be defined!")
@test_throws err BBMaterial(; horizon = 1, rho = 8e-6, E = 210e3)

mat4 = BBMaterial(horizon=1.5, rho=7850.0, E=210e9, Gc=1.0)
io = IOBuffer()
show(io, "text/plain", mat4)
msg = String(take!(io))

correct_msg1 = "BBMaterial:\n  δ:   1.5\n  rho: 7850.0\n  E:   2.1e11\n  "
correct_msg1 *= "nu:  0.25\n  G:   8.4e10\n  K:   1.4e11\n  bc:  1.584475877892647e11\n  "
correct_msg1 *= "Gc:  1.0\n  εc:  1.6265001215808886e-6"

correct_msg2 = "Peridynamics.BBMaterial:\n  δ:   1.5\n  rho: 7850.0\n  E:   2.1e11\n  "
correct_msg2 *= "nu:  0.25\n  G:   8.4e10\n  K:   1.4e11\n  bc:  1.584475877892647e11\n  "
correct_msg2 *= "Gc:  1.0\n  εc:  1.6265001215808886e-6"

@test msg == correct_msg1 || msg == correct_msg2
