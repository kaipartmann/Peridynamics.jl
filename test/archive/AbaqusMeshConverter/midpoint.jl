using Test
using Peridynamics.AbaqusMeshConverter: midpoint

@test midpoint([1], [2], [3], [4]) == [2.5]
@test midpoint(([0, 0, 0] for _ in 1:4)...) == [0, 0, 0]
@test midpoint(([0, 0, 0] for _ in 1:8)...) == [0, 0, 0]
@test midpoint(([i, i+1, i+2] for i in 1:4)...) == [2.5, 3.5, 4.5]
@test midpoint(([i, i+1, i+2] for i in 1:8)...) == [4.5, 5.5, 6.5]
@test_throws MethodError midpoint([1], [2], [3])
