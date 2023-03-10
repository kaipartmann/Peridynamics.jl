using Test
using Peridynamics.AbaqusMeshConverter: tetvol

@test tetvol([0,0,0], [1,0,0], [0,1,0], [0,0,1]) == 1/6
@test tetvol(([0,0,0] for _ in 1:4)...) == 0
