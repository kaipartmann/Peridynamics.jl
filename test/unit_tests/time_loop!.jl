using Peridynamics, Test

if Threads.nthreads() <= 2
    positions = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    point_spacing = 1.0
    δ = 1.5 * point_spacing
    n_points = 2
    volumes = fill(point_spacing^3, n_points)
    pc = PointCloud(positions, volumes)
    mat = BondBasedMaterial(; horizon=δ, rho=1, E=1, Gc=1)
    body = Peridynamics.create_simmodel(mat, pc)
    sim = PDSingleBodyAnalysis(;
        name="",
        pc=pc,
        mat=mat,
        bcs=[VelocityBC(t -> 0.1, [2], 1)],
        td=VelocityVerlet(1, 1),
        es=ExportSettings(),
    )
    Peridynamics.time_loop!(body, sim.td, sim.mat, sim.bcs, sim.ics, sim.es)

    @test body.velocity_half == [0.0 0.1; 0.0 0.0; 0.0 0.0]
    @test body.displacement == [0.0 0.1; 0.0 0.0; 0.0 0.0]
    @test body.position == [0.0 1.1; 0.0 0.0; 0.0 0.0]

    b¹² = [
        18 * mat.E / (3 * (1 - 2 * 0.25)) / (π * δ^4) * 1.1 * 0.1/1.1 * 1.0
        0.0
        0.0
    ]
    b_int = [b¹² -b¹²]
    @test body.b_int[:,:,1] ≈ b_int
    @test body.acceleration ≈ b_int
    @test body.velocity ≈ [0.0 0.1; 0.0 0.0; 0.0 0.0] + b_int * 0.5

else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end

if Threads.nthreads() <= 2
    positions = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    point_spacing = 1.0
    δ = 1.5 * point_spacing
    n_points = 2
    volumes = fill(point_spacing^3, n_points)
    pc = PointCloud(positions, volumes)
    mat = BondBasedMaterial(; horizon=δ, rho=1, E=1, Gc=1)
    body = Peridynamics.create_simmodel(mat, pc)
    sim = PDSingleBodyAnalysis(;
        name="",
        pc=pc,
        mat=mat,
        bcs=[ForceDensityBC(t -> -1, [1], 1), ForceDensityBC(t -> 1, [2], 1)],
        td=DynamicRelaxation(1),
        es=ExportSettings(),
    )
    Peridynamics.time_loop!(body, sim.td, sim.mat, sim.bcs, sim.ics, sim.es)

    @test body.velocity_half == [-0.09375 0.09375; 0.0 0.0; 0.0 0.0]
    @test body.displacement == zeros(3, 2)
    @test body.position == positions
    @test body.b_int == zeros(3, 2, Threads.nthreads())
    @test body.acceleration == zeros(3, 2)
    @test body.velocity == [-0.046875 0.046875; 0.0 0.0; 0.0 0.0]
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
