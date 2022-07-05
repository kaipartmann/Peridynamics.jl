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
    @test body.n_active_family_members == [1; 1;;]
    @test body.bond_failure == [1]

    # Boundary Condition:
    # Point 2 with ẋ_z = 1 m/s with Δt = 0.0015 s
    body.position[1, 2] = 1.0015
    Peridynamics.compute_forcedensity!(body, mat)

    @test body.n_active_family_members == [0; 0;;]
    @test body.bond_failure == [0]

    b¹² = [
        18 * E / (3 * (1 - 2 * 0.25)) / (π * δ^4) * 1.0015 * 0.0015/1.0015 * 1.0
        0.0
        0.0
    ]
    Threads.@threads for _ in 1:Threads.nthreads()
        for i in body.owned_points[Threads.threadid()]
            body.b_int[1,i,1] = sum(@view body.b_int[1,i,:])
            body.b_int[2,i,1] = sum(@view body.b_int[2,i,:])
            body.b_int[3,i,1] = sum(@view body.b_int[3,i,:])
        end
    end
    @test body.b_int[:,:,1] ≈ [b¹² -b¹²]
else
    @warn "Test omitted! Threads.nthreads() should be <= 2"
end
