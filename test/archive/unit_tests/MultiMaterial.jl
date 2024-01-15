using Peridynamics, Test

mat1 = BBMaterial(horizon = 1, rho = 1, E = 1, Gc = 1)
mat2 = BBMaterial(horizon = 2, rho = 2, E = 2, Gc = 2)
materials = (mat1, mat2)
pointmap = [1, 2, 1, 2]
mm = Peridynamics.MultiMaterial(materials, pointmap)

for p in pointmap
    @test mm[p] == materials[p]
end
@test mm[:] == (mat1, mat2, mat1, mat2)
@test mm[begin:end] == (mat1, mat2, mat1, mat2)
@test mm[[2, 3]] == (mat2, mat1)
@test eltype(mm) == BBMaterial

io = IOBuffer()
show(io, "text/plain", mm)
msg = String(take!(io))
@test msg == "MultiMaterial{BBMaterial, 2, UInt8} for 4 points" ||
      msg ==
      "Peridynamics.MultiMaterial{Peridynamics.BBMaterial, 2, UInt8} for 4 points"

err_msg = "N = 1\nUsage of MultiMaterial only makes sense for N > 1"
@test_throws ErrorException(err_msg) Peridynamics.MultiMaterial((mat1,), [1, 2, 3, 4, 5, 6])

err_msg = "Index in `matofpoint` out of bounds!"
@test_throws ErrorException(err_msg) Peridynamics.MultiMaterial(materials, [1, 2, 3])

# many materials -> I != UInt8
N = typemax(UInt8) + 1
materials = Tuple(mat1 for _ in 1:N)
pointmap = rand(1:N, N)
mm = Peridynamics.MultiMaterial(materials, pointmap)
@test typeof(mm) == MultiMaterial{BBMaterial, 256, Int64}

pointmap = convert(Vector{UInt16}, pointmap)
mm = Peridynamics.MultiMaterial(materials, pointmap)
@test typeof(mm) == MultiMaterial{BBMaterial, 256, UInt16}

# Interface for AbstractMaterial:
@test mat1[1] == mat1
@test mat1[1:3] == mat1
@test mat1[:] == mat1
@test firstindex(mat1) == 1
@test lastindex(mat1) == 1
@test eltype(mat1) == BBMaterial
