@testitem "bond-associated family" begin
    # setup
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BACCMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, nu=0.25, Gc=1)
    pd = Peridynamics.PointDecomposition(body, 1)

    # 1
    system = Peridynamics.BondAssociatedSystem(body, pd, 1)

    @test system.position == position
    @test system.volume == volume
    @test system.bonds == [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(3, √2, true),
        Peridynamics.Bond(4, √2, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(2, √2, true),
        Peridynamics.Bond(4, √2, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(2, √2, true),
        Peridynamics.Bond(3, √2, true),
    ]
    @test system.n_neighbors == [3, 3, 3, 3]
    @test system.bond_ids == [1:3, 4:6, 7:9, 10:12]

end

@testitem "intersection bond ids" begin
    pos, vol = uniform_box(1, 0.25, 0.25, 0.25)
    body = Body(BACCMaterial(), pos, vol)
    material!(body, horizon=0.5, rho=1, E=1, nu=0.25, Gc=1)

    pd = Peridynamics.PointDecomposition(body, 1)

    bas = Peridynamics.get_system(body, pd, 1)
    (; bonds, bond_ids, intersection_bond_ids) = bas

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

    for i in Peridynamics.each_point_idx(system)
        for bond_idx in Peridynamics.each_bond_idx(system, i)
            # bond = bonds[bond_idx]
            # j = bond.neighbor
            for babond_idx in Peridynamics.each_intersecting_bond_idx(system, i, bond_idx)

            end
        end
    end

end
