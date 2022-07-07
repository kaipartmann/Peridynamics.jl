using Peridynamics
using Test

##------------------------------------------------------------------------------------------
# stable time step calculation

if Threads.nthreads() <= 2
    positions = [
        0.0 1.0
        0.0 0.0
        0.0 0.0
    ]
    point_spacing = 1.0
    δ = 1.5 * point_spacing
    n_points = 2
    volumes = fill(point_spacing^3, n_points)
    pc = PointCloud(positions, volumes)
    mat = BondBasedMaterial(; horizon=δ, rho=1, E=1, Gc=1)
    body = Peridynamics.create_simmodel(mat, pc)

    Δt = 0.7 * sqrt(2 * 1 / (body.volumes[2] * 1 / 1 * 18 * 2/3 / (π * δ^4)))
    @test Peridynamics.calc_stable_timestep(body, mat.rho, mat.K, δ) == Δt
    @test calc_stable_user_timestep(pc, mat) == Δt
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end

##------------------------------------------------------------------------------------------
# Velocity verlet

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
        td=TimeDiscretization(1, 1),
        es=ExportSettings(),
    )
    Peridynamics.velocity_verlet!(body, sim)

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


##------------------------------------------------------------------------------------------
# Adaptive dynamic relaxation

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
        td=TimeDiscretization(1; alg=:dynrelax),
        es=ExportSettings(),
    )
    Peridynamics.dynamic_relaxation_finite_difference!(body, sim)

    @test body.velocity_half == [-0.0625 0.0625; 0.0 0.0; 0.0 0.0]
    @test body.displacement == [-0.0625 0.0625; 0.0 0.0; 0.0 0.0]
    @test body.position == [-0.0625 1.0625; 0.0 0.0; 0.0 0.0]
    @test body.b_int == zeros(3, 2, 1)
    @test body.acceleration == zeros(3, 2)
    @test body.velocity == [-0.03125 0.03125; 0.0 0.0; 0.0 0.0]
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
