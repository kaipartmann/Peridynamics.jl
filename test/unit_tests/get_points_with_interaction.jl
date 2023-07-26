using Peridynamics, Test

if Threads.nthreads() == 2
    pc = PointCloud(1, 0.5, 0.5, 0.25)
    pc.failure_flag .= false
    mat1 = ContinuumBasedMaterial(; horizon=0.8, rho=8000, E=210e9, nu=0.2, Gc=1, C1=1, C2=1, C3=1)
    mat2 = deepcopy(mat1)
    pointmap = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2]
    mat = MultiMaterial((mat1, mat2), pointmap)
    owned_points = Peridynamics.defaultdist(pc.n_points, Threads.nthreads())

    owned_points_with_twoni = Peridynamics.points_with_twoni(mat, owned_points)
    @test owned_points_with_twoni == owned_points

    owned_points_with_threeni = Peridynamics.points_with_threeni(mat, owned_points)
    @test owned_points_with_threeni == owned_points

    mat2 = ContinuumBasedMaterial(; horizon=0.8, rho=8000, E=210e9, nu=0.2, Gc=1, C1=1, C2=0, C3=0)
    mat = MultiMaterial((mat1, mat2), pointmap)
    owned_points = Peridynamics.defaultdist(pc.n_points, Threads.nthreads())

    owned_points_with_twoni = Peridynamics.points_with_twoni(mat, owned_points)
    @test owned_points_with_twoni == [[1,2,3,4],[5,6,7,8]]

    owned_points_with_threeni = Peridynamics.points_with_threeni(mat, owned_points)
    @test owned_points_with_threeni == [[1,2,3,4],[5,6,7,8]]
else
    @warn "Test omitted! Threads.nthreads() should be == 2"
end