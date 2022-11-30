using Peridynamics
using Peridynamics: defaultdist, find_bonds, find_unique_bonds, create_simmodel,
    define_precrack!, calc_damage!
using Test

##------------------------------------------------------------------------------------------
# PointCloud
@testset "PointCloud" begin
    position1 = [
        1.0 2.0
        3.0 4.0
        5.0 6.0
    ]
    volume1 = [1.0, 1.0]
    pc1 = PointCloud(position1, volume1)
    @test pc1.n_points == 2
    @test pc1.failure_flag == [true, true]
    @test pc1.radius == Peridynamics.sphere_radius.(volume1)
    @test isempty(pc1.point_sets)
    @test_throws DimensionMismatch PointCloud([1.0 2.0; 3.0 4.0], volume1)

    pc2 = PointCloud(1, 1, 1, 0.5)
    @test pc2.position == [
        -0.25  -0.25  -0.25  -0.25   0.25   0.25   0.25  0.25
        -0.25  -0.25   0.25   0.25  -0.25  -0.25   0.25  0.25
        -0.25   0.25  -0.25   0.25  -0.25   0.25  -0.25  0.25
    ]
    @test pc2.n_points == 8
    @test pc2.volume == fill(0.125, 8)
    @test pc2.radius == Peridynamics.sphere_radius.(fill(0.125, 8))
    @test pc2.failure_flag == fill(true, 8)
    @test isempty(pc2.point_sets)

    pc3 = PointCloud(1, 1, 1, 0.5; center_x = 1, center_y = 1, center_z = 1)
    @test pc3.position == [
        0.75  0.75  0.75  0.75  1.25  1.25  1.25  1.25
        0.75  0.75  1.25  1.25  0.75  0.75  1.25  1.25
        0.75  1.25  0.75  1.25  0.75  1.25  0.75  1.25
    ]
    @test pc3.n_points == 8
    @test pc3.volume == fill(0.125, 8)
    @test pc3.radius == Peridynamics.sphere_radius.(fill(0.125, 8))
    @test pc3.failure_flag == fill(true, 8)
    @test isempty(pc3.point_sets)
end

##------------------------------------------------------------------------------------------
# bonds

@testset "bonds" begin
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
        owned_points = defaultdist(n_points, Threads.nthreads())
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
        body = create_simmodel(mat, pc)
        precrack = PreCrack([1], [2])

        if Threads.nthreads() == 1
            @test body.n_active_family_members == [1; 1;;]
        elseif Threads.nthreads() == 2
            @test body.n_active_family_members == [1 0; 0 1]
        end
        @test body.bond_failure == [1]
        @test body.damage == [0, 0]

        define_precrack!(body, precrack)
        calc_damage!(body)

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
end

@testset "show PointCloud" begin
    positions = [
        0.0 1.0
        0.0 0.0
        0.0 0.0
    ]
    volumes = ones(2)
    pc = PointCloud(positions, volumes)
    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (20, 40)), "text/plain", pc)
    msg_pc = String(take!(io))
    @test msg_pc == "2-points PointCloud"
end
