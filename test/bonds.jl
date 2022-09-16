using Peridynamics
using Peridynamics: defaultdist, find_bonds, find_unique_bonds
using Test

if Threads.nthreads() <= 2
    positions = [
        0.0 1.0
        0.0 0.0
        0.0 0.0
    ]
    n_points = 2
    point_spacing = 1.0
    δ = 1.5 * point_spacing
    pc = PointCloud(positions, ones(n_points))
    owned_points = Peridynamics.defaultdist(n_points, Threads.nthreads())
    bond_data, n_family_members = find_bonds(pc, δ, owned_points)

    @test bond_data == [(1,2,1.,true),(2,1,1.,true)]
    @test n_family_members == [1,1]

    bond_data, n_family_members = find_bonds(pc, 0.9, owned_points)

    @test bond_data == Vector{Tuple{Int,Int,Float64,Bool}}()
    @test n_family_members == [0,0]

    bond_data, n_family_members = find_unique_bonds(pc, 1.1, owned_points)

    @test bond_data == [(1,2,1.,true)]
    @test n_family_members == [1,1]

    bond_data, n_family_members = find_unique_bonds(pc,0.9,owned_points)

    @test bond_data == Vector{Tuple{Int,Int,Float64,Bool}}()
    @test n_family_members == [0,0]

    mat = BondBasedMaterial(horizon=δ, rho=1, E=1, Gc=1)
    body = Peridynamics.create_simmodel(mat, pc)
    precrack = PreCrack([1], [2])

    if Threads.nthreads() == 1
        @test body.n_active_family_members == [1; 1;;]
    elseif Threads.nthreads() == 2
        @test body.n_active_family_members == [1 0; 0 1]
    end
    @test body.bond_failure == [1]
    @test body.damage == [0, 0]

    Peridynamics.define_precrack!(body, precrack)
    Peridynamics.calc_damage!(body)

    if Threads.nthreads() == 1
        @test body.n_active_family_members == [0; 0;;]
    elseif Threads.nthreads() == 2
        @test body.n_active_family_members == [0 0; 0 0]
    end
    @test body.bond_failure == [0]
    @test body.damage == [1, 1]
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
