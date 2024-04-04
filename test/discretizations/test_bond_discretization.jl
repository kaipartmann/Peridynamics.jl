@testitem "init_bond_system" begin
    # setup
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    pd = Peridynamics.PointDecomposition(body, 2)

    # 1
    bd, ch = Peridynamics.init_bond_system(body, pd, 1)

    @test bd.position == position
    @test bd.volume == volume
    @test bd.bonds == [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(3, √2, true),
        Peridynamics.Bond(4, √2, true),
    ]
    @test bd.n_neighbors == [3, 3]
    @test bd.bond_ids == [1:3, 4:6]

    @test ch.point_ids == [1, 2, 3, 4]
    @test ch.loc_points == [1, 2]
    @test ch.halo_points == [3, 4]
    @test ch.hidxs_by_src[2] == 3:4

    for i in 1:4
        @test ch.localizer[i] == i
    end

    # 2
    bd, ch = Peridynamics.init_bond_system(body, pd, 2)

    @test bd.position == position[:, [3, 4, 1, 2]]
    @test bd.volume == volume[[3, 4, 1, 2]]
    @test bd.bonds == [
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, √2, true),
        Peridynamics.Bond(2, √2, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, √2, true),
        Peridynamics.Bond(1, √2, true),
    ]
    @test bd.n_neighbors == [3, 3]
    @test bd.bond_ids == [1:3, 4:6]

    @test ch.point_ids == [3, 4, 1, 2]
    @test ch.loc_points == [3, 4]
    @test ch.halo_points == [1, 2]
    @test ch.hidxs_by_src[1] == 3:4
    @test ch.localizer[3] == 1
    @test ch.localizer[4] == 2
    @test ch.localizer[1] == 3
    @test ch.localizer[2] == 4
end

@testitem "find_bonds!" begin
    # setup
    position = [0.0 1.0
                0.0 0.0
                0.0 0.0]
    fail_permit = [true, true]
    balltree = Peridynamics.NearestNeighbors.BallTree(position)

    # find point 2
    δ = 1.5
    bonds = Vector{Peridynamics.Bond}()
    n_neighbors = Peridynamics.find_bonds!(bonds, balltree, position, fail_permit, δ, 1)
    @test n_neighbors == 1
    @test bonds == [Peridynamics.Bond(2, 1.0, true)]

    # horizon too small - find nothing
    δ = 0.9
    bonds = Vector{Peridynamics.Bond}()
    n_neighbors = Peridynamics.find_bonds!(bonds, balltree, position, fail_permit, δ, 1)
    @test n_neighbors == 0
    @test bonds == Vector{Peridynamics.Bond}()

    # no failure allowed for point 2
    fail_permit[2] = false
    δ = 1.5
    bonds = Vector{Peridynamics.Bond}()
    n_neighbors = Peridynamics.find_bonds!(bonds, balltree, position, fail_permit, δ, 1)
    @test n_neighbors == 1
    @test bonds == [Peridynamics.Bond(2, 1.0, false)]
end

@testitem "find_bonds" begin
    # setup
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1, 1, 1, 1]
    mat = BBMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)

    # all points are local points
    loc_points = 1:4
    bonds, n_neighbors = Peridynamics.find_bonds(body, loc_points)
    @test bonds == [
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
    @test n_neighbors == [3, 3, 3, 3]

    loc_points = 1:2
    bonds, n_neighbors = Peridynamics.find_bonds(body, loc_points)
    @test bonds == [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(3, √2, true),
        Peridynamics.Bond(4, √2, true),
    ]
    @test n_neighbors == [3, 3]

    loc_points = 2:3
    bonds, n_neighbors = Peridynamics.find_bonds(body, loc_points)
    @test bonds == [
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(3, √2, true),
        Peridynamics.Bond(4, √2, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(2, √2, true),
        Peridynamics.Bond(4, √2, true),
    ]
    @test n_neighbors == [3, 3]
end

@testitem "find_halo_points" begin
    bonds = [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
    ]

    # no halo point
    loc_points = 1:4
    halo_points = Peridynamics.find_halo_points(bonds, loc_points)
    @test halo_points == Int[]

    # only 1 halo point
    loc_points = 1:3
    halo_points = Peridynamics.find_halo_points(bonds, loc_points)
    @test halo_points == [4]

    # 2 halo points
    loc_points = 1:2
    halo_points = Peridynamics.find_halo_points(bonds, loc_points)
    @test halo_points == [3, 4] || halo_points == [4, 3]

    # 2 halo points
    loc_points = 3:4
    halo_points = Peridynamics.find_halo_points(bonds, loc_points)
    @test halo_points == [1, 2] || halo_points == [2, 1]
end

@testitem "find_bond_ids" begin
    n_neighbors = [1, 2]
    bond_ids = Peridynamics.find_bond_ids(n_neighbors)
    @test bond_ids == [1:1, 2:3]

    n_neighbors = [3, 4, 5]
    bond_ids = Peridynamics.find_bond_ids(n_neighbors)
    @test bond_ids == [1:3, 4:7, 8:12]
end
