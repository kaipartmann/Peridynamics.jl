using Peridynamics, Test

if Threads.nthreads() <= 2
    positions = [
        0.0 1.0
        0.0 0.0
        0.0 0.0
    ]
    n_points = 2
    pc = PointCloud(positions, ones(n_points))
    mat = BondBasedMaterial(horizon=1.5, rho=1, E=1, Gc=1)
    mat_δ⁻ = BondBasedMaterial(horizon=0.9, rho=1, E=1, Gc=1)

    owned_points = Peridynamics.defaultdist(n_points, Threads.nthreads())
    bond_data, n_family_members = Peridynamics.find_bonds(pc, mat, owned_points)

    @test bond_data == [(1,2,1.,true),(2,1,1.,true)]
    @test n_family_members == [1,1]

    bond_data, n_family_members = Peridynamics.find_bonds(pc, mat_δ⁻, owned_points)

    @test bond_data == Vector{Tuple{Int,Int,Float64,Bool}}()
    @test n_family_members == [0,0]

    bond_data, n_family_members = Peridynamics.find_unique_bonds(pc, mat, owned_points)

    @test bond_data == [(1,2,1.,true)]
    @test n_family_members == [1,1]

    bond_data, n_family_members = Peridynamics.find_unique_bonds(pc, mat_δ⁻, owned_points)

    @test bond_data == Vector{Tuple{Int,Int,Float64,Bool}}()
    @test n_family_members == [0,0]

    body = Peridynamics.init_body(mat, pc)
    precrack = PreCrack([1], [2])

    if Threads.nthreads() == 1
        @test body.n_active_family_members == [1; 1;;]
    elseif Threads.nthreads() == 2
        @test body.n_active_family_members == [1 0; 1 0]
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
