using Peridynamics, Test

let
    if Threads.nthreads() <= 3
        # Test one- and two-neighbour-interactions, 3 points
        position = [0.0 1.0 0.0
                    0.0 0.0 0.0
                    0.0 0.0 1.0]
        volume = fill(1.0, 3)
        pc = PointCloud(position, volume)
        mat = CPDMaterial(; horizon=1.5, rho=8e-6, E=210e3, nu=0.3, Gc=1.0, C1=1, C2=1, C3=1)
        owned_points = Peridynamics.defaultdist(pc.n_points, Threads.nthreads())
        bond_data, n_family_members = Peridynamics.find_bonds(pc, mat, owned_points)
        hood_range = Peridynamics.find_hood_range(n_family_members, pc.n_points)
        two_ni_data, n_two_ni, n_two_ni_per_point = Peridynamics.find_two_ni(pc, mat, bond_data, hood_range, owned_points)

        @test bond_data == [(1, 2, 1.0, true),
                            (1, 3, 1.0, true),
                            (2, 1, 1.0, true),
                            (2, 3, √2, true),
                            (3, 1, 1.0, true),
                            (3, 2, √2, true)]
        @test n_family_members == [2,2,2]
        @test hood_range == [1:2,3:4,5:6]
        @test two_ni_data == [(1,2,1.0),
                              (3,4,1.0),
                              (5,6,1.0)]
        @test n_two_ni == 3
        @test n_two_ni_per_point == [1,1,1]
    else
        @warn "Test omitted! Threads.nthreads() should be <= 3"
    end
end

let
    if Threads.nthreads() <= 5
        # Test two- and three-neigbour-interactions, 5 points
        position = [0.0 1.0 0.0 0.0 2.0
                    0.0 0.0 1.0 0.0 2.0
                    0.0 0.0 0.0 1.0 2.0]
        volume = fill(1.0, 5)
        pc = PointCloud(position, volume)
        mat = CPDMaterial(; horizon=1.5, rho=8e-6, E=210e3, nu=0.3, Gc=1.0, C1=1, C2=1, C3=1)
        owned_points = Peridynamics.defaultdist(pc.n_points, Threads.nthreads())
        bond_data, n_family_members = Peridynamics.find_bonds(pc, mat, owned_points)
        hood_range = Peridynamics.find_hood_range(n_family_members, pc.n_points)
        two_ni_data, n_two_ni, n_two_ni_per_point = Peridynamics.find_two_ni(pc, mat, bond_data, hood_range, owned_points)
        three_ni_data, n_three_ni, n_three_ni_per_point = Peridynamics.find_three_ni(pc, mat, bond_data, hood_range, owned_points)

        @test bond_data == [(1, 2, 1.0, true),
                            (1, 3, 1.0, true),
                            (1, 4, 1.0, true),
                            (2, 1, 1.0, true),
                            (2, 3, √2, true),
                            (2, 4, √2, true),
                            (3, 1, 1.0, true),
                            (3, 2, √2, true),
                            (3, 4, √2, true),
                            (4, 1, 1.0, true),
                            (4, 2, √2, true),
                            (4, 3, √2, true)]
        @test n_family_members == [3,3,3,3,0]
        @test hood_range == [1:3,4:6,7:9,10:12,13:12]
        @test two_ni_data == [(1,2,1.0),
                              (1,3,1.0),
                              (2,3,1.0),
                              (4,5,1.0),
                              (4,6,1.0),
                              (5,6,√3),
                              (7,8,1.0),
                              (7,9,1.0),
                              (8,9,√3),
                              (10,11,1.0),
                              (10,12,1.0),
                              (11,12,√3)]
        @test n_two_ni == 12
        @test n_two_ni_per_point == [3,3,3,3,0]
        @test three_ni_data == [(1,2,3,1.0),
                                (4,5,6,1.0),
                                (7,8,9,1.0),
                                (10,11,12,1.0)]
        @test n_three_ni == 4
        @test n_three_ni_per_point == [1,1,1,1,0]
    else
        @warn "Test omitted! Threads.nthreads() should be <= 5"
    end
end
