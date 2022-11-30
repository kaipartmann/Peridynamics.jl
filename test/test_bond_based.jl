using Peridynamics
using Test

##------------------------------------------------------------------------------------------
# Bond based material

mat1 = BondBasedMaterial(; horizon = 1, rho = 1, E = 1, Gc = 1)
@test mat1.δ == 1
@test mat1.rho == 1
@test mat1.E == 1
@test mat1.nu == 1/4
@test mat1.G == 1 / (2 * (1 + 1/4))
@test mat1.K ≈ 1 / (3 * (1 - 2 * 1/4))
@test mat1.bc ≈ 18 * 1 / (3 * (1 - 2 * 1/4)) / (π * 1^4)
@test mat1.Gc == 1
@test mat1.εc == sqrt(5.0 * 1 / (9.0 * 1 / (3 * (1 - 2 * 1/4)) * 1))

mat2 = BondBasedMaterial(; horizon = 0, rho = 0, E = 0, Gc = 0)
@test mat2.δ == 0
@test mat2.rho == 0
@test mat2.E == 0
@test mat2.nu == 1/4
@test mat2.G == 0
@test mat2.K == 0
@test isnan(mat2.bc)
@test mat2.Gc == 0
@test isnan(mat2.εc)

##------------------------------------------------------------------------------------------
# Force density calculation and BondBasedBody creation

if Threads.nthreads() <= 2
    positions = [
        0.0 1.0
        0.0 0.0
        0.0 0.0
    ]
    point_spacing = 1.0
    δ = 1.5 * point_spacing
    E = 210e9
    n_points = 2
    volumes = fill(point_spacing^3, n_points)
    pc = PointCloud(positions, volumes)
    mat = BondBasedMaterial(horizon=δ, rho=7850.0, E=E, Gc=1.0)
    body = Peridynamics.create_simmodel(mat, pc)

    @test body.n_points == n_points
    @test body.n_bonds == 1
    @test body.n_threads == Threads.nthreads()
    @test body.unique_bonds == true
    @test body.b_int == zeros(3, n_points, Threads.nthreads())
    if Threads.nthreads() == 1
        @test body.n_active_family_members == [1; 1;;]
    elseif Threads.nthreads() == 2
        @test body.n_active_family_members == [1 0; 0 1]
    end
    @test body.bond_failure == [1]

    # Boundary Condition:
    # Point 2 with ẋ_z = 1 m/s with Δt = 0.0015 s
    body.position[1, 2] = 1.0015
    Peridynamics.compute_forcedensity!(body, mat)

    if Threads.nthreads() == 1
        @test body.n_active_family_members == [0; 0;;]
    elseif Threads.nthreads() == 2
        @test body.n_active_family_members == [0 0; 0 0]
    end
    @test body.bond_failure == [0]

    b¹² = [
        18 * E / (3 * (1 - 2 * 0.25)) / (π * δ^4) * 1.0015 * 0.0015/1.0015 * 1.0
        0.0
        0.0
    ]
    Peridynamics.update_thread_cache!(body)
    @test body.b_int[:,:,1] ≈ [b¹² -b¹²]
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end

if Threads.nthreads() <= 2
    positions = [
        0.0 1.0
        0.0 0.0
        0.0 0.0
    ]
    point_spacing = 1.0
    δ = 1.5 * point_spacing
    E = 210e9
    n_points = 2
    volumes = fill(point_spacing^3, n_points)
    pc = PointCloud(positions, volumes)
    mat = BondBasedMaterial(horizon=δ, rho=7850.0, E=E, Gc=1.0)
    body = Peridynamics.create_simmodel(mat, pc)

    # show BondBasedMaterial
    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (20, 40)), "text/plain", mat)
    msg_mat = String(take!(io))
    @test msg_mat == string(
        "BondBasedMaterial:\n  δ:   1.5\n  rho: 7850.0\n  E:   2.1e11\n  ",
        "nu:  0.25\n  G:   8.4e10\n  K:   1.4e11\n  bc:  1.584475877892647e11\n  ",
        "Gc:  1.0\n  εc:  1.6265001215808886e-6",
    )

    # show BondBasedBody
    io = IOBuffer()
    show(IOContext(io, :limit => true, :displaysize => (20, 40)), "text/plain", body)
    msg_body = String(take!(io))
    @test msg_body == "2-points Peridynamics.BondBasedBody:\n  1 bonds"
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
