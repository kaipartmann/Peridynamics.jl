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
    pc = PointCloud(positions, ones(n_points))
    owned_points = Peridynamics.defaultdist(n_points, Threads.nthreads())
    one_ni_data, n_family_members = find_bonds(pc, 1.1, owned_points)

    @test one_ni_data == [(1,2,1.,true),(2,1,1.,true)]
    @test n_family_members == [1,1]

    one_ni_data, n_family_members = find_bonds(pc,0.9,owned_points)

    @test one_ni_data == Vector{Tuple{Int,Int,Float64,Bool}}()
    @test n_family_members == [0,0]

    one_ni_data, n_family_members = find_unique_bonds(pc, 1.1, owned_points)

    @test one_ni_data == [(1,2,1.,true)]
    @test n_family_members == [1,1]

    one_ni_data, n_family_members = find_unique_bonds(pc,0.9,owned_points)

    @test one_ni_data == Vector{Tuple{Int,Int,Float64,Bool}}()
    @test n_family_members == [0,0]
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
