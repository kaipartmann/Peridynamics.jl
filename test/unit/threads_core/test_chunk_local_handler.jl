using Peridynamics, Test

## point cloud with only two points
let
    # setup
    position = [
        0.0 1.0
        0.0 0.0
        0.0 0.0
    ]
    n_points = 2
    pc = Peridynamics.PointCloud(position, ones(n_points))
    mat = Peridynamics.BBMaterial(horizon=1.5, rho=1, E=1, Gc=1)

    # n_chunks = 1
    vclh = Peridynamics.chunk_local_handlers(pc, mat, 1)
    @test length(vclh) == 1
    clh = first(vclh)
    @test clh.point_ids == [1, 2]
    @test clh.loc_points == 1:2
    @test clh.halo_points == Int[]
    @test clh.localizer == Dict(2=>2, 1=>1)
    @test typeof(clh.d) == Peridynamics.PointBondDiscretization
    @test clh.d.position == position
    @test clh.d.volume == ones(n_points)
    @test clh.d.bonds == [Peridynamics.Bond(2, 1.0, true), Peridynamics.Bond(1,1.0,true)]
    @test clh.d.n_neighbors == [1, 1]
    @test clh.d.bond_range == [1:1, 2:2]
    @test typeof(clh.s) == Peridynamics.BBStorage
    @test clh.s.position == position
    @test clh.s.displacement == zeros(3, n_points)
    @test clh.s.velocity == zeros(3, n_points)
    @test clh.s.velocity_half == zeros(3, n_points)
    @test clh.s.acceleration == zeros(3, n_points)
    @test clh.s.b_int == zeros(3, n_points)
    @test clh.s.b_ext == zeros(3, n_points)
    @test clh.s.damage == zeros(n_points)
    @test clh.s.bond_active == [true, true]
    @test clh.s.n_active_bonds == [1, 1]

    # n_chunks = 2
    vclh = Peridynamics.chunk_local_handlers(pc, mat, 2)
    @test length(vclh) == 2
    clh = vclh[1]
    @test clh.point_ids == [1, 2]
    @test clh.loc_points == 1:1
    @test clh.halo_points == [2]
    @test clh.localizer == Dict(2=>2, 1=>1)
    @test typeof(clh.d) == Peridynamics.PointBondDiscretization
    @test clh.d.position == position
    @test clh.d.volume == ones(n_points)
    @test clh.d.bonds == [Peridynamics.Bond(2, 1.0, true)]
    @test clh.d.n_neighbors == [1]
    @test clh.d.bond_range == [1:1]
    @test typeof(clh.s) == Peridynamics.BBStorage
    @test clh.s.position == position
    @test clh.s.displacement == zeros(3, 1)
    @test clh.s.velocity == zeros(3, 1)
    @test clh.s.velocity_half == zeros(3, 1)
    @test clh.s.acceleration == zeros(3, 1)
    @test clh.s.b_int == zeros(3, 1)
    @test clh.s.b_ext == zeros(3, 1)
    @test clh.s.damage == zeros(1)
    @test clh.s.bond_active == [true]
    @test clh.s.n_active_bonds == [1]
    clh = vclh[2]
    @test clh.point_ids == [2, 1]
    @test clh.loc_points == 2:2
    @test clh.halo_points == [1]
    @test clh.localizer == Dict(2=>1, 1=>2)
    @test typeof(clh.d) == Peridynamics.PointBondDiscretization
    @test clh.d.position == [1.0 0.0; 0.0 0.0; 0.0 0.0]
    @test clh.d.volume == ones(n_points)
    @test clh.d.bonds == [Peridynamics.Bond(2, 1.0, true)] # index is localized!
    @test clh.d.n_neighbors == [1]
    @test clh.d.bond_range == [1:1]
    @test typeof(clh.s) == Peridynamics.BBStorage
    @test clh.s.position == [1.0 0.0; 0.0 0.0; 0.0 0.0]
    @test clh.s.displacement == zeros(3, 1)
    @test clh.s.velocity == zeros(3, 1)
    @test clh.s.velocity_half == zeros(3, 1)
    @test clh.s.acceleration == zeros(3, 1)
    @test clh.s.b_int == zeros(3, 1)
    @test clh.s.b_ext == zeros(3, 1)
    @test clh.s.damage == zeros(1)
    @test clh.s.bond_active == [true]
    @test clh.s.n_active_bonds == [1]
end

## point cloud with 4 points
let
    # setup
    position = [
        0.0 1.0 2.0 3.0
        0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0
    ]
    n_points = 4
    pc = Peridynamics.PointCloud(position, ones(n_points))
    mat = Peridynamics.BBMaterial(horizon=1.1, rho=1, E=1, Gc=1)

    # n_chunks = 2
    vclh = Peridynamics.chunk_local_handlers(pc, mat, 2)
    @test length(vclh) == 2
    clh = vclh[1]
    @test clh.point_ids == [1, 2, 3]
    @test clh.loc_points == 1:2
    @test clh.halo_points == [3]
    @test clh.localizer == Dict(1=>1, 2=>2, 3=>3)
    @test typeof(clh.d) == Peridynamics.PointBondDiscretization
    @test clh.d.position == [0.0 1.0 2.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    @test clh.d.volume == [1.0, 1.0, 1.0]
    @test clh.d.bonds == [Peridynamics.Bond(2, 1.0, true), Peridynamics.Bond(1, 1.0, true),
                          Peridynamics.Bond(3, 1.0, true)]
    @test clh.d.n_neighbors == [1,2]
    @test clh.d.bond_range == [1:1, 2:3]
    @test typeof(clh.s) == Peridynamics.BBStorage
    @test clh.s.position == [0.0 1.0 2.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    @test clh.s.displacement == zeros(3, 2)
    @test clh.s.velocity == zeros(3, 2)
    @test clh.s.velocity_half == zeros(3, 2)
    @test clh.s.acceleration == zeros(3, 2)
    @test clh.s.b_int == zeros(3, 2)
    @test clh.s.b_ext == zeros(3, 2)
    @test clh.s.damage == zeros(2)
    @test clh.s.bond_active == [true, true, true]
    @test clh.s.n_active_bonds == [1, 2]
    clh = vclh[2]
    @test clh.point_ids == [3, 4, 2]
    @test clh.loc_points == 3:4
    @test clh.halo_points == [2]
    @test clh.localizer == Dict(3=>1, 4=>2, 2=>3)
    @test typeof(clh.d) == Peridynamics.PointBondDiscretization
    @test clh.d.position == [2.0 3.0 1.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    @test clh.d.volume == [1.0, 1.0, 1.0]
    @test clh.d.bonds == [Peridynamics.Bond(3, 1.0, true), Peridynamics.Bond(2, 1.0, true),
                          Peridynamics.Bond(1, 1.0, true)]
    @test clh.d.n_neighbors == [2, 1]
    @test clh.d.bond_range == [1:2, 3:3]
    @test typeof(clh.s) == Peridynamics.BBStorage
    @test clh.s.position == [2.0 3.0 1.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    @test clh.s.displacement == zeros(3, 2)
    @test clh.s.velocity == zeros(3, 2)
    @test clh.s.velocity_half == zeros(3, 2)
    @test clh.s.acceleration == zeros(3, 2)
    @test clh.s.b_int == zeros(3, 2)
    @test clh.s.b_ext == zeros(3, 2)
    @test clh.s.damage == zeros(2)
    @test clh.s.bond_active == [true, true, true]
    @test clh.s.n_active_bonds == [2, 1]
end
