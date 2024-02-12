using Peridynamics, Test
import Peridynamics: PointCloud

## find_bonds!
let
    # setup
    position = [
        0.0 1.0
        0.0 0.0
        0.0 0.0
    ]
    n_points = 2
    pc = PointCloud(position, ones(n_points))

    # find point 2
    bond_buf = Vector{Peridynamics.Bond}()
    n_neighbors = Peridynamics.find_bonds!(bond_buf, pc, 1.5, 1)
    @test n_neighbors == 1
    @test bond_buf == [Peridynamics.Bond(2, 1.0, true)]

    # horizon too small - find nothing
    bond_buf = Vector{Peridynamics.Bond}()
    n_neighbors = Peridynamics.find_bonds!(bond_buf, pc, 0.9, 1)
    @test n_neighbors == 0
    @test bond_buf == Vector{Peridynamics.Bond}()

    # no failure allowed for point 2
    pc.failure_allowed[2] = false
    bond_buf = Vector{Peridynamics.Bond}()
    n_neighbors = Peridynamics.find_bonds!(bond_buf, pc, 1.5, 1)
    @test n_neighbors == 1
    @test bond_buf == [Peridynamics.Bond(2, 1.0, false)]
end

## find_bonds with BBMaterial
let
    # setup
    position = [
        0.0 1.0 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 0.0 0.0 1.0
    ]
    n_points = 4
    pc = Peridynamics.PointCloud(position, ones(n_points))
    mat = Peridynamics.BBMaterial(horizon=1.5, rho=1, E=1, Gc=1)

    # all points are local points
    loc_points = 1:4
    bonds, n_neighbors = Peridynamics.find_bonds(pc, mat, loc_points)
    @test bonds == [
        Peridynamics.Bond(2, 1.0, true)
        Peridynamics.Bond(3, 1.0, true)
        Peridynamics.Bond(4, 1.0, true)
        Peridynamics.Bond(1, 1.0, true)
        Peridynamics.Bond(3, √2, true)
        Peridynamics.Bond(4, √2, true)
        Peridynamics.Bond(1, 1.0, true)
        Peridynamics.Bond(2, √2, true)
        Peridynamics.Bond(4, √2, true)
        Peridynamics.Bond(1, 1.0, true)
        Peridynamics.Bond(2, √2, true)
        Peridynamics.Bond(3, √2, true)
    ]
    @test n_neighbors == [3, 3, 3, 3]
end

## localize!
let
    # change two bonds
    bonds = [Peridynamics.Bond(100, 1.0, true), Peridynamics.Bond(101, 30.0, false)]
    localizer = Dict(100 => 1, 101 => 2)
    Peridynamics.localize!(bonds, localizer)
    @test bonds == [Peridynamics.Bond(1, 1.0, true), Peridynamics.Bond(2, 30.0, false)]

    # do not change any bond
    bonds = [Peridynamics.Bond(100, 1.0, true), Peridynamics.Bond(101, 30.0, false)]
    localizer = Dict(100 => 100, 101 => 101)
    Peridynamics.localize!(bonds, localizer)
    @test bonds == [Peridynamics.Bond(100, 1.0, true), Peridynamics.Bond(101, 30.0, false)]

    # key not found error
    bonds = [Peridynamics.Bond(100, 1.0, true)]
    localizer = Dict(2 => 1)
    @test_throws KeyError(100) Peridynamics.localize!(bonds, localizer)
end

## find_bond_range
let
    n_neighbors = [1, 2]
    bond_range = Peridynamics.find_bond_range(n_neighbors)
    @test bond_range == [1:1, 2:3]

    n_neighbors = [3, 4, 5]
    bond_range = Peridynamics.find_bond_range(n_neighbors)
    @test bond_range == [1:3, 4:7, 8:12]
end

## find_halo_points
let
    bonds = [
        Peridynamics.Bond(2, 1.0, true)
        Peridynamics.Bond(3, 1.0, true)
        Peridynamics.Bond(4, 1.0, true)
        Peridynamics.Bond(1, 1.0, true)
    ]

    # only 1 halo point
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
    @test halo_points == [3,4] || halo_points == [4,3]

    # 2 halo points
    loc_points = 3:4
    halo_points = Peridynamics.find_halo_points(bonds, loc_points)
    @test halo_points == [1,2] || halo_points == [2,1]
end
