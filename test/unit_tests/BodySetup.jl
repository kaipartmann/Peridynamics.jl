using Peridynamics, Test

bs1 = BodySetup(PointCloud(1, 1, 1, 0.5), BondBasedMaterial(; horizon=1, rho=1, E=1, Gc=1))
@test isempty(bs1.precracks)
@test isempty(bs1.bcs)
@test isempty(bs1.ics)
