using Peridynamics, Test

mat1 = BondBasedMaterial(horizon=1, rho=1, E=1, Gc=1)
mat2 = BondBasedMaterial(horizon=2, rho=2, E=2, Gc=2)
materials = (mat1, mat2)
pointmap = [1,2,1,2]
mm = Peridynamics.MultiMaterial(materials, pointmap)

for p in pointmap
    @test mm[p] == materials[p]
end
@test mm[:] == (mat1, mat2, mat1, mat2)
@test mm[begin:end] == (mat1, mat2, mat1, mat2)
@test mm[[2, 3]] == (mat2, mat1)
@test eltype(mm) == BondBasedMaterial

io = IOBuffer()
show(io, "text/plain", mm)
msg = String(take!(io))
@test msg == "MultiMaterial{BondBasedMaterial, 2} for 4 points" ||
      msg == "Peridynamics.MultiMaterial{Peridynamics.BondBasedMaterial, 2} for 4 points"
