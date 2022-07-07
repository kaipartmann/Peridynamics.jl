using Peridynamics
using Test

##------------------------------------------------------------------------------------------
# defaultdist

@test Peridynamics.defaultdist(4, 2) == [1:2, 3:4]
@test Peridynamics.defaultdist(3, 2) == [1:2, 3:3]
@test Peridynamics.defaultdist(2, 3) == [1:1, 2:2, 0:0]
