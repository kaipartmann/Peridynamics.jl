@testitem "intersection bond ids" begin
    pos, vol = uniform_box(1, 0.25, 0.25, 0.25)
    body = Body(BANOSBMaterial(), pos, vol)
    material!(body, horizon=0.5, rho=1, E=1, nu=0.25, Gc=1)

    pd = Peridynamics.PointDecomposition(body, 1)

    bas, ch = Peridynamics.get_system(body, pd, 1)
    (; bonds, bond_ids, intersection_bond_ids, volume_factor) = bas

    @test bonds == [
        Peridynamics.Bond(2, 0.25, true)
        Peridynamics.Bond(3, 0.5, true)
        Peridynamics.Bond(1, 0.25, true)
        Peridynamics.Bond(3, 0.25, true)
        Peridynamics.Bond(4, 0.5, true)
        Peridynamics.Bond(1, 0.5, true)
        Peridynamics.Bond(2, 0.25, true)
        Peridynamics.Bond(4, 0.25, true)
        Peridynamics.Bond(2, 0.5, true)
        Peridynamics.Bond(3, 0.25, true)
    ]
    @test bond_ids == [1:2, 3:5, 6:8, 9:10]
    @test intersection_bond_ids == [[2], [1], [], [3], [2], [2], [1], [], [2], [1]]
    @test volume_factor ≈ [
        0.5
        0.5
        0.0
        0.3333333333333333
        0.3333333333333333
        0.3333333333333333
        0.3333333333333333
        0.0
        0.5
        0.5
    ]
end
