@testitem "DOF handling, raw functions" begin
    A = zeros(Int, 3, 10)
    for (dof, dim, i) in Peridynamics.each_dof_idx(size(A, 1), axes(A, 2))
        A[dof] = dof
    end
    @test A[1, 1] == 1
    @test A[2, 1] == 2
    @test A[3, 1] == 3
    @test A[1, 2] == 4
    @test A[2, 2] == 5
    @test A[3, 2] == 6
    @test A[1, 10] == 28
    @test A[2, 10] == 29
    @test A[3, 10] == 30

    @test Peridynamics.get_dof(3, 1, 1) == 1
    @test Peridynamics.get_dof(3, 2, 1) == 2
    @test Peridynamics.get_dof(3, 3, 1) == 3
    @test Peridynamics.get_dof(3, 1, 2) == 4
    @test Peridynamics.get_dof(3, 2, 2) == 5
    @test Peridynamics.get_dof(3, 3, 2) == 6
    @test Peridynamics.get_dof(3, 1, 10) == 28
    @test Peridynamics.get_dof(3, 2, 10) == 29
    @test Peridynamics.get_dof(3, 3, 10) == 30

    @test Peridynamics.get_point(3, 1) == 1
    @test Peridynamics.get_point(3, 2) == 1
    @test Peridynamics.get_point(3, 3) == 1
    @test Peridynamics.get_point(3, 4) == 2
    @test Peridynamics.get_point(3, 5) == 2
    @test Peridynamics.get_point(3, 6) == 2
    @test Peridynamics.get_point(3, 28) == 10
    @test Peridynamics.get_point(3, 29) == 10
    @test Peridynamics.get_point(3, 30) == 10

    @test Peridynamics.get_dim(3, 1) == 1
    @test Peridynamics.get_dim(3, 2) == 2
    @test Peridynamics.get_dim(3, 3) == 3
    @test Peridynamics.get_dim(3, 4) == 1
    @test Peridynamics.get_dim(3, 5) == 2
    @test Peridynamics.get_dim(3, 6) == 3
    @test Peridynamics.get_dim(3, 28) == 1
    @test Peridynamics.get_dim(3, 29) == 2
    @test Peridynamics.get_dim(3, 30) == 3
end

@testitem "DOF handling, BondSystem interface, -t 1" begin
    position, volume = uniform_box(1,1,1,0.25)
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    ts = VelocityVerlet(steps=10)
    pd = Peridynamics.PointDecomposition(body, 1)
    system = Peridynamics.get_system(body, pd, 1)

    @test Peridynamics.get_n_dim(system) == 3
    @test Peridynamics.get_n_points(system) == 64
    @test Peridynamics.get_n_loc_points(system) == 64
    @test Peridynamics.get_n_dof(system) == 192
    @test Peridynamics.get_n_loc_dof(system) == 192
    @test Peridynamics.get_dof(system, 1, 1) == 1
    @test Peridynamics.get_dof(system, 2, 1) == 2
    @test Peridynamics.get_dof(system, 3, 1) == 3
    @test Peridynamics.get_dof(system, 1, 2) == 4
    @test Peridynamics.get_dof(system, 2, 2) == 5
    @test Peridynamics.get_dof(system, 3, 2) == 6
    @test Peridynamics.get_dof(system, 1, 10) == 28
    @test Peridynamics.get_dof(system, 2, 10) == 29
    @test Peridynamics.get_dof(system, 3, 10) == 30
    @test Peridynamics.get_dof(system, 1, 64) == 190

    @test Peridynamics.each_dim(system) == 1:3
    all_dof_idxs = collect(Peridynamics.each_dof_idx(system))
    @test size(all_dof_idxs) == (64, 3)
    @test all_dof_idxs[1, 1] == (1, 1, 1)
    @test all_dof_idxs[1, 2] == (2, 2, 1)
    @test all_dof_idxs[1, 3] == (3, 3, 1)
    @test all_dof_idxs[2, 1] == (4, 1, 2)
    @test all_dof_idxs[2, 2] == (5, 2, 2)
    @test all_dof_idxs[2, 3] == (6, 3, 2)
    @test all_dof_idxs[10, 1] == (28, 1, 10)
    @test all_dof_idxs[10, 2] == (29, 2, 10)
    @test all_dof_idxs[10, 3] == (30, 3, 10)
    @test all_dof_idxs[64, 1] == (190, 1, 64)
    @test all_dof_idxs == collect(Peridynamics.each_loc_dof_idx(system))
    @test collect(Peridynamics.each_dof_idx(system, [1,2,10])) == [
        (1, 1, 1) (2, 2, 1) (3, 3, 1)
        (4, 1, 2) (5, 2, 2) (6, 3, 2)
        (28, 1, 10) (29, 2, 10) (30, 3, 10)
    ]
    all_dofs = collect(Peridynamics.each_dof(system))
    @test size(all_dofs) == (64, 3)
    @test all_dofs[1, 1] == 1
    @test all_dofs[1, 2] == 2
    @test all_dofs[1, 3] == 3
    @test all_dofs[2, 1] == 4
    @test all_dofs[2, 2] == 5
    @test all_dofs[2, 3] == 6
    @test all_dofs[10, 1] == 28
    @test all_dofs[10, 2] == 29
    @test all_dofs[10, 3] == 30
    @test all_dofs[64, 1] == 190
    @test all_dofs == collect(Peridynamics.each_loc_dof(system))

    @test Peridynamics.get_point(system, 1) == 1
    @test Peridynamics.get_dim(system, 1) == 1
    @test Peridynamics.get_point(system, 2) == 1
    @test Peridynamics.get_dim(system, 2) == 2
    @test Peridynamics.get_point(system, 3) == 1
    @test Peridynamics.get_dim(system, 3) == 3
    @test Peridynamics.get_point(system, 4) == 2
    @test Peridynamics.get_dim(system, 4) == 1
    @test Peridynamics.get_point(system, 28) == 10
    @test Peridynamics.get_dim(system, 28) == 1
    @test Peridynamics.get_point(system, 29) == 10
    @test Peridynamics.get_dim(system, 29) == 2
    @test Peridynamics.get_point(system, 30) == 10
    @test Peridynamics.get_dim(system, 30) == 3
end

@testitem "DOF handling, BondSystem interface, -t 2" begin
    position, volume = uniform_box(1,1,1,0.5)
    body = Body(BBMaterial(), position, volume)
    material!(body, horizon=2, rho=1, E=1, Gc=1)
    ts = VelocityVerlet(steps=10)
    pd = Peridynamics.PointDecomposition(body, 2)

    X1 = [-0.25, -0.25, -0.25]
    X4 = [0.25, 0.25, -0.25]
    X5 = [-0.25, -0.25, 0.25]
    X8 = [0.25, 0.25, 0.25]

    # first chunk
    system = Peridynamics.get_system(body, pd, 1)

    # point ids on the first chunk similar to the body point ids
    @test system.position[:, 1] ≈ X1
    @test system.position[:, 4] ≈ X4
    @test system.position[:, 5] ≈ X5
    @test system.position[:, 8] ≈ X8

    @test Peridynamics.get_n_dim(system) == 3
    @test Peridynamics.get_n_points(system) == 8
    @test Peridynamics.get_n_loc_points(system) == 4
    @test Peridynamics.get_n_dof(system) == 24
    @test Peridynamics.get_n_loc_dof(system) == 12
    @test Peridynamics.get_dof(system, 1, 1) == 1
    @test Peridynamics.get_dof(system, 2, 1) == 2
    @test Peridynamics.get_dof(system, 3, 1) == 3
    @test Peridynamics.get_dof(system, 1, 2) == 4
    @test Peridynamics.get_dof(system, 2, 2) == 5
    @test Peridynamics.get_dof(system, 3, 2) == 6
    @test Peridynamics.get_dof(system, 1, 8) == 22
    @test Peridynamics.get_dof(system, 2, 8) == 23
    @test Peridynamics.get_dof(system, 3, 8) == 24

    @test Peridynamics.each_dim(system) == 1:3
    all_dof_idxs = collect(Peridynamics.each_dof_idx(system))
    @test size(all_dof_idxs) == (8, 3)
    @test all_dof_idxs[1, 1] == (1, 1, 1)
    @test all_dof_idxs[1, 2] == (2, 2, 1)
    @test all_dof_idxs[1, 3] == (3, 3, 1)
    @test all_dof_idxs[2, 1] == (4, 1, 2)
    @test all_dof_idxs[2, 2] == (5, 2, 2)
    @test all_dof_idxs[2, 3] == (6, 3, 2)
    @test all_dof_idxs[8, 1] == (22, 1, 8)
    @test all_dof_idxs[8, 2] == (23, 2, 8)
    @test all_dof_idxs[8, 3] == (24, 3, 8)
    @test all_dof_idxs[1:4, :] == collect(Peridynamics.each_loc_dof_idx(system))
    @test collect(Peridynamics.each_dof_idx(system, [1,2,8])) == [
        (1, 1, 1) (2, 2, 1) (3, 3, 1)
        (4, 1, 2) (5, 2, 2) (6, 3, 2)
        (22, 1, 8) (23, 2, 8) (24, 3, 8)
    ]
    all_dofs = collect(Peridynamics.each_dof(system))
    @test size(all_dofs) == (8, 3)
    @test all_dofs[1, 1] == 1
    @test all_dofs[1, 2] == 2
    @test all_dofs[1, 3] == 3
    @test all_dofs[2, 1] == 4
    @test all_dofs[2, 2] == 5
    @test all_dofs[2, 3] == 6
    @test all_dofs[8, 1] == 22
    @test all_dofs[8, 2] == 23
    @test all_dofs[8, 3] == 24
    @test all_dofs[1:4, :] == collect(Peridynamics.each_loc_dof(system))

    @test Peridynamics.get_point(system, 1) == 1
    @test Peridynamics.get_dim(system, 1) == 1
    @test Peridynamics.get_point(system, 2) == 1
    @test Peridynamics.get_dim(system, 2) == 2
    @test Peridynamics.get_point(system, 3) == 1
    @test Peridynamics.get_dim(system, 3) == 3
    @test Peridynamics.get_point(system, 4) == 2
    @test Peridynamics.get_dim(system, 4) == 1
    # no bounds checking, so the following still works although we have only 8 points
    @test Peridynamics.get_point(system, 28) == 10
    @test Peridynamics.get_dim(system, 28) == 1
    @test Peridynamics.get_point(system, 29) == 10
    @test Peridynamics.get_dim(system, 29) == 2
    @test Peridynamics.get_point(system, 30) == 10
    @test Peridynamics.get_dim(system, 30) == 3

    # second chunk
    system = Peridynamics.get_system(body, pd, 2)

    # point ids on the second chunk NOT similar to the body point ids
    @test system.position[:, 1] ≈ X5
    @test system.position[:, 4] ≈ X8
    @test system.position[:, 5] ≈ X1
    @test system.position[:, 8] ≈ X4

    @test Peridynamics.get_n_dim(system) == 3
    @test Peridynamics.get_n_points(system) == 8
    @test Peridynamics.get_n_loc_points(system) == 4
    @test Peridynamics.get_n_dof(system) == 24
    @test Peridynamics.get_n_loc_dof(system) == 12
    @test Peridynamics.get_dof(system, 1, 1) == 1
    @test Peridynamics.get_dof(system, 2, 1) == 2
    @test Peridynamics.get_dof(system, 3, 1) == 3
    @test Peridynamics.get_dof(system, 1, 2) == 4
    @test Peridynamics.get_dof(system, 2, 2) == 5
    @test Peridynamics.get_dof(system, 3, 2) == 6
    @test Peridynamics.get_dof(system, 1, 8) == 22
    @test Peridynamics.get_dof(system, 2, 8) == 23
    @test Peridynamics.get_dof(system, 3, 8) == 24

    @test Peridynamics.each_dim(system) == 1:3
    all_dof_idxs = collect(Peridynamics.each_dof_idx(system))
    @test size(all_dof_idxs) == (8, 3)
    @test all_dof_idxs[1, 1] == (1, 1, 1)
    @test all_dof_idxs[1, 2] == (2, 2, 1)
    @test all_dof_idxs[1, 3] == (3, 3, 1)
    @test all_dof_idxs[2, 1] == (4, 1, 2)
    @test all_dof_idxs[2, 2] == (5, 2, 2)
    @test all_dof_idxs[2, 3] == (6, 3, 2)
    @test all_dof_idxs[8, 1] == (22, 1, 8)
    @test all_dof_idxs[8, 2] == (23, 2, 8)
    @test all_dof_idxs[8, 3] == (24, 3, 8)
    @test all_dof_idxs[1:4, :] == collect(Peridynamics.each_loc_dof_idx(system))
    @test collect(Peridynamics.each_dof_idx(system, [1,2,8])) == [
        (1, 1, 1) (2, 2, 1) (3, 3, 1)
        (4, 1, 2) (5, 2, 2) (6, 3, 2)
        (22, 1, 8) (23, 2, 8) (24, 3, 8)
    ]
    all_dofs = collect(Peridynamics.each_dof(system))
    @test size(all_dofs) == (8, 3)
    @test all_dofs[1, 1] == 1
    @test all_dofs[1, 2] == 2
    @test all_dofs[1, 3] == 3
    @test all_dofs[2, 1] == 4
    @test all_dofs[2, 2] == 5
    @test all_dofs[2, 3] == 6
    @test all_dofs[8, 1] == 22
    @test all_dofs[8, 2] == 23
    @test all_dofs[8, 3] == 24
    @test all_dofs[1:4, :] == collect(Peridynamics.each_loc_dof(system))

    @test Peridynamics.get_point(system, 1) == 1
    @test Peridynamics.get_dim(system, 1) == 1
    @test Peridynamics.get_point(system, 2) == 1
    @test Peridynamics.get_dim(system, 2) == 2
    @test Peridynamics.get_point(system, 3) == 1
    @test Peridynamics.get_dim(system, 3) == 3
    @test Peridynamics.get_point(system, 4) == 2
    @test Peridynamics.get_dim(system, 4) == 1
    # no bounds checking, so the following still works although we have only 8 points
    @test Peridynamics.get_point(system, 28) == 10
    @test Peridynamics.get_dim(system, 28) == 1
    @test Peridynamics.get_point(system, 29) == 10
    @test Peridynamics.get_dim(system, 29) == 2
    @test Peridynamics.get_point(system, 30) == 10
    @test Peridynamics.get_dim(system, 30) == 3
end
