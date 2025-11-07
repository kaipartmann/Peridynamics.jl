@testitem "damage changed flag" begin
    pos, vol = uniform_box(1, 1, 1, 0.4)
    body = Body(RKCMaterial(), pos, vol)
    material!(body; horizon=1.5, rho=1, E=210e9, nu=0.25, Gc=1.0)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; n_neighbors) = chunk.system
    (; update_gradients, damage, n_active_bonds, bond_active, gradient_weight) = chunk.storage

    @test n_neighbors == fill(7, 8)
    @test n_active_bonds == fill(7, 8)
    @test bond_active == fill(true, 56)
    @test gradient_weight[:,1] ≈ [1.1548446667165417, -0.3316495750857463, -0.3316495750857462]
    @test gradient_weight[:,2] ≈ [-0.33164957508574605, 1.1548446667165417, -0.3316495750857461]
    @test update_gradients == fill(false, 8)

    Peridynamics.calc_weights_and_defgrad!(chunk, 0.0, 0.0)

    # everything should be the same, the gradients are initialized and damage did not change
    @test n_active_bonds == fill(7, 8)
    @test bond_active == fill(true, 56)
    @test gradient_weight[:,1] ≈ [1.1548446667165417, -0.3316495750857463, -0.3316495750857462]
    @test gradient_weight[:,2] ≈ [-0.33164957508574605, 1.1548446667165417, -0.3316495750857461]
    @test update_gradients == fill(false, 8)

    # bond 1 failed somehow
    bond_active[1] = false
    Peridynamics.calc_weights_and_defgrad!(chunk, 0.0, 0.0)

    # now the changes should be reflected in the chunk
    @test n_active_bonds == [6, 7, 7, 7, 7, 7, 7, 7]
    @test bond_active == [false; fill(true, 55)]
    @test iszero(gradient_weight[:,1])
    @test gradient_weight[:,2] ≈ [-0.6163778391975847, 1.2366132461014017, -0.24988099570088682]
    @test update_gradients == fill(false, 8)
end

@testitem "RKCMaterial initialization" begin
    # Test default constructor
    mat1 = RKCMaterial()
    @test mat1.kernel == cubic_b_spline_kernel_norm
    @test mat1.constitutive_model isa SaintVenantKirchhoff
    @test mat1.dmgmodel isa CriticalStretch
    @test mat1.monomial == :C1
    @test mat1.lambda == 0
    @test mat1.beta ≈ sqrt(eps())

    # Test constructor with parameters
    mat2 = RKCMaterial(
        kernel = linear_kernel,
        model = LinearElastic(),
        dmgmodel = CriticalStretch(),
        monomial = :C1,
        lambda = 0,
        beta = 1e-8
    )
    @test mat2.kernel == linear_kernel
    @test mat2.constitutive_model isa LinearElastic
    @test mat2.dmgmodel isa CriticalStretch
    @test mat2.monomial == :C1
    @test mat2.lambda == 0
    @test mat2.beta == 1e-8

    # Test constructor with invalid lambda/beta
    @test_throws ArgumentError RKCMaterial(lambda = 2.0)
    @test_throws ArgumentError RKCMaterial(lambda = -0.5)
    @test_throws ArgumentError RKCMaterial(beta = 2.0)
    @test_throws ArgumentError RKCMaterial(beta = -0.5)
end

@testitem "RKCRMaterial initialization" begin
    # Test default constructor
    mat1 = RKCRMaterial()
    @test mat1.kernel == cubic_b_spline_kernel_norm
    @test mat1.constitutive_model isa SaintVenantKirchhoff
    @test mat1.dmgmodel isa CriticalStretch
    @test mat1.monomial == :C1
    @test mat1.lambda == 0
    @test mat1.beta ≈ sqrt(eps())

    # Test constructor with parameters (only linear elastic models are supported)
    mat2 = RKCRMaterial(
        kernel = linear_kernel,
        model = LinearElastic(),
        dmgmodel = CriticalStretch(),
        monomial = :C1,
        lambda = 0,
        beta = 1e-10,
    )
    @test mat2.kernel == linear_kernel
    @test mat2.constitutive_model isa LinearElastic
    @test mat2.dmgmodel isa CriticalStretch
    @test mat2.monomial == :C1
    @test mat2.lambda == 0
    @test mat2.beta == 1e-10

    # Test failure with non-LinearElastic model
    @test_throws ArgumentError RKCRMaterial(model = NeoHooke())

    # Test constructor with invalid lambda/beta
    @test_throws ArgumentError RKCRMaterial(lambda = 2.0)
    @test_throws ArgumentError RKCRMaterial(lambda = -0.5)
    @test_throws ArgumentError RKCRMaterial(beta = 2.0)
    @test_throws ArgumentError RKCRMaterial(beta = -0.5)
end

@testitem "gradient weights calculation" begin
    using Peridynamics.LinearAlgebra

    pos, vol = uniform_box(1, 1, 1, 0.4)
    body = Body(RKCMaterial(), pos, vol)
    material!(body; horizon=1.5, rho=1, E=210e9, nu=0.25, Gc=1.0)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; gradient_weight, weighted_volume, position) = chunk.storage

    # Reset and recalculate weights
    gradient_weight .= 0.0
    for i in 1:8
        Peridynamics.rkc_weights!(chunk.storage, chunk.system, chunk.mat, chunk.paramsetup, 0.0, 0.0, i)
    end

    # Test that gradient weights are consistent and correctly calculated
    @test gradient_weight[:,1] ≈ [1.1548446667165417, -0.3316495750857463, -0.3316495750857462]
    @test gradient_weight[:,2] ≈ [-0.33164957508574605, 1.1548446667165417, -0.3316495750857461]
    @test gradient_weight[:,3] ≈ [0.5612621576497658, 0.5612621576497656, -0.45224360054794666]
    @test gradient_weight[:,4] ≈ [-0.3316495750857459, -0.3316495750857462, 1.1548446667165417]
    @test gradient_weight[:,5] ≈ [0.5612621576497658, -0.45224360054794643, 0.561262157649766]
    @test gradient_weight[:,6] ≈ [-0.45224360054794627, 0.5612621576497654, 0.5612621576497657]
    @test gradient_weight[:,7] ≈ [0.22263101798392687, 0.2226310179839266, 0.22263101798392676]
    @test gradient_weight[:,8] ≈ [-1.1548446667165413, -0.33164957508574583, -0.331649575085746]
    @test weighted_volume[1] > 0
end

@testitem "deformation gradient calculation for RKCMaterial" begin
    using Peridynamics.StaticArrays, Peridynamics.LinearAlgebra

    pos, vol = uniform_box(1, 1, 1, 0.1)
    body = Body(RKCMaterial(), pos, vol)
    material!(body; horizon=1.5, rho=1, E=210e9, nu=0.25, Gc=1.0)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; storage, system, paramsetup, mat) = chunk
    (; position, displacement, defgrad) = storage

    # No displacement, should be identity matrix
    Peridynamics.calc_force_density!(dh, 0.0, 0.0)

    for i in 1:64
        F = Peridynamics.get_tensor(defgrad, i)
        @test F ≈ [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] atol=1e-12
    end

    # Apply small uniform stretch in x-direction
    F_a = @SMatrix [1.00001 0.0 0.0
                    0.0 1.0 0.0
                    0.0 0.0 1.0]
    for i in eachindex(vol)
        position[:,i] = F_a * position[:,i]
    end

    Peridynamics.calc_force_density!(dh, 0.0, 0.0)

    # Check that deformation gradient captures the stretch
    for i in eachindex(vol)
        F = Peridynamics.get_tensor(defgrad, i)
        @test F ≈ F_a atol=1e-5
    end
end

@testitem "deformation gradient calculation for RKCRMaterial" begin
    using Peridynamics.StaticArrays, Peridynamics.LinearAlgebra

    pos, vol = uniform_box(1, 1, 1, 0.1)
    body = Body(RKCRMaterial(), pos, vol)
    material!(body; horizon=1.5, rho=1, E=210e9, nu=0.25, Gc=1.0)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; storage, system, paramsetup, mat) = chunk
    (; position, displacement, defgrad, defgrad_dot, velocity_half) = storage

    # No displacement, should be identity matrix
    Peridynamics.calc_force_density!(dh, 0.0, 0.0)

    for i in eachindex(vol)
        F = Peridynamics.get_tensor(defgrad, i)
        @test F ≈ [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] atol=1e-12
    end

    # Apply small uniform stretch in x-direction
    F_a = @SMatrix [1.00001 0.0 0.0
                    0.0 1.0 0.0
                    0.0 0.0 1.0]
    for i in eachindex(vol)
        X = Peridynamics.get_vector(position, i)
        X_new = F_a * X
        Peridynamics.update_vector!(position, i, X_new)
    end

    Peridynamics.calc_force_density!(dh, 0.0, 0.0)

    # Check that deformation gradient captures the stretch
    for i in eachindex(vol)
        F = Peridynamics.get_tensor(defgrad, i)
        @test F ≈ F_a atol=1e-5
    end

    # Add a velocity field to check defgrad_dot
    velocity_half .= 0.0
    for i in eachindex(vol)
        velocity_half[1,i] = position[1,i] * 0.1
    end

    # Recalculate with velocity
    Peridynamics.calc_force_density!(dh, 0.0, 0.0)

    # Check that defgrad_dot is calculated correctly
    for i in eachindex(vol)
        F_dot = Peridynamics.get_tensor(defgrad_dot, i)
        @test F_dot[1,1] > 0  # Should have positive rate in x direction
        @test abs(F_dot[2,2]) < 1e-14
        @test abs(F_dot[3,3]) < 1e-14
    end
end

# @testitem "stress calculation RKCMaterial" begin
#     using Peridynamics.StaticArrays, Peridynamics.LinearAlgebra

#     pos, vol = uniform_box(1, 1, 1, 0.1)
#     body = Body(RKCMaterial(), pos, vol)
#     material!(body; horizon=1.5, rho=1, E=210e9, nu=0.25, Gc=1.0)

#     dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
#     chunk = dh.chunks[1]
#     (; storage, system, paramsetup, mat) = chunk
#     # (; position, displacement, defgrad, defgrad_dot, velocity_half, stress) = storage
#     (; position, displacement, defgrad, velocity_half, stress) = storage
#     (; first_piola_kirchhoff) = storage

#     # Apply small uniform stretch in x-direction
#     F_a = @SMatrix [1.00001 0.0 0.0
#                     0.0 1.0 0.0
#                     0.0 0.0 1.0]
#     for i in eachindex(vol)
#         X = Peridynamics.get_vector(position, i)
#         X_new = F_a * X
#         Peridynamics.update_vector!(position, i, X_new)
#     end

#     # Calculate force density (which calculates stress)
#     Peridynamics.calc_force_density!(dh, 0.0, 0.0)

#     # Check stress tensor
#     for i in 1:8
#         σ = Peridynamics.get_tensor(first_piola_kirchhoff, i)
#         # For small strain, stress should primarily be in xx component
#         @test σ[1,1] > 0  # Positive stress in x direction
#         @test abs(σ[2,3]) < 1e-10  # No shear in yz
#         @test abs(σ[1,3]) < 1e-10  # No shear in xz
#         @test abs(σ[1,2]) < 1e-10  # No shear in xy
#     end
# end

# @testitem "stress calculation RKCRMaterial with rotation" begin
#     using Peridynamics.LinearAlgebra

#     pos, vol = uniform_box(1, 1, 1, 0.4)
#     body = Body(RKCRMaterial(), pos, vol)
#     material!(body; horizon=1.5, rho=1, E=210e9, nu=0.25, Gc=1.0)

#     dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
#     chunk = dh.chunks[1]
#     (; defgrad, position, displacement, rotation, left_stretch, bond_stress) = chunk.storage

#     # Apply both stretch and rotation
#     stretch_factor = 1.001
#     theta = π/60  # Small rotation angle

#     for i in 1:8
#         x, y = position[1,i], position[2,i]
#         # Apply rotation + stretch
#         x_new = (cos(theta) * x - sin(theta) * y) * stretch_factor
#         y_new = (sin(theta) * x + cos(theta) * y)

#         displacement[:,i] = [x_new - x, y_new - y, 0.0]
#         position[:,i] += displacement[:,i]
#     end

#     # Calculate force density (which calculates stress)
#     Peridynamics.calc_force_density!(dh, 0.0, 0.0)

#     # Check rotation matrix and left stretch for at least one bond
#     bond_id = 1
#     R = Peridynamics.get_tensor(rotation, bond_id)
#     V = Peridynamics.get_tensor(left_stretch, bond_id)

#     # The determinant of R should be close to 1
#     @test abs(det(R) - 1.0) < 1e-10

#     # R should be orthogonal (R * R' = I)
#     @test R * transpose(R) ≈ I atol=1e-10

#     # V should be symmetric
#     @test V ≈ transpose(V) atol=1e-10

#     # Check bond stress
#     σ = Peridynamics.get_tensor(bond_stress, bond_id)
#     @test σ[1,1] > 0  # Positive stress in x due to stretch
# end

@testitem "damage model integration with RKCMaterial" begin
    pos, vol = uniform_box(1, 1, 1, 0.4)
    body = Body(RKCMaterial(dmgmodel=CriticalStretch()), pos, vol)
    material!(body; horizon=1.5, rho=1, E=210e9, nu=0.25, Gc=1.0)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; bond_active, damage, n_active_bonds, position, displacement) = chunk.storage

    # Apply large stretch to exceed critical stretch
    stretch_factor = 1.01  # Should be large enough to break bonds

    # Apply stretch to points on one side only to create a strain
    for i in 5:8  # Assuming points 5-8 are on one side
        displacement[:,i] = [position[1,i] * (stretch_factor - 1), 0.0, 0.0]
        position[:,i] += displacement[:,i]
    end

    # Calculate force density with damage evaluation
    Peridynamics.calc_force_density!(dh, 0.0, 0.0)

    # Verify some bonds are broken
    @test count(bond_active) < length(bond_active)

    # Verify damage values are calculated
    @test any(damage .> 0.0)

    # Verify n_active_bonds reflects broken bonds
    @test any(n_active_bonds .< 7)
end

@testitem "reproducing kernel basis functions" begin
    using Peridynamics.StaticArrays, Peridynamics.LinearAlgebra

    # Test C1 kernel
    @test Peridynamics.get_q_dim(:C1) == 3
    ΔX = [0.1, 0.2, 0.3]
    Q_C1 = Peridynamics.get_monomial_vector(:C1, ΔX)
    @test Q_C1 == SVector(0.1, 0.2, 0.3)
    Q∇ᵀ_C1 = Peridynamics.get_gradient_extraction_matrix(:C1)
    @test Q∇ᵀ_C1 == SMatrix{3,3}([1 0 0; 0 1 0; 0 0 1])

    # Test RK1 kernel
    @test Peridynamics.get_q_dim(:RK1) == 4
    Q_RK1 = Peridynamics.get_monomial_vector(:RK1, ΔX)
    @test Q_RK1 == SVector(1.0, 0.1, 0.2, 0.3)
    Q∇ᵀ_RK1 = Peridynamics.get_gradient_extraction_matrix(:RK1)
    @test Q∇ᵀ_RK1 == SMatrix{3,4}([0 1 0 0; 0 0 1 0; 0 0 0 1])

    # Test RK2 kernel
    @test Peridynamics.get_q_dim(:RK2) == 7
    Q_RK2 = Peridynamics.get_monomial_vector(:RK2, ΔX)
    @test Q_RK2 ≈ SVector(1.0, 0.1, 0.2, 0.3, 0.01, 0.04, 0.09) atol=1e-15
    Q∇ᵀ_RK2 = Peridynamics.get_gradient_extraction_matrix(:RK2)
    expected_RK2 = SMatrix{3,7}([0 1 0 0 0 0 0
                                 0 0 1 0 0 0 0;
                                 0 0 0 1 0 0 0])
    @test Q∇ᵀ_RK2 == expected_RK2

    # Test PD2 kernel
    @test Peridynamics.get_q_dim(:PD2) == 9
    Q_PD2 = Peridynamics.get_monomial_vector(:PD2, ΔX)
    @test Q_PD2 ≈ SVector(0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 0.04, 0.06, 0.09) atol=1e-15
    Q∇ᵀ_PD2 = Peridynamics.get_gradient_extraction_matrix(:PD2)
    expected_PD2 = SMatrix{3,9}([1 0 0 0 0 0 0 0 0;
                                 0 1 0 0 0 0 0 0 0;
                                 0 0 1 0 0 0 0 0 0])
    @test Q∇ᵀ_PD2 == expected_PD2

    # Test error handling for unknown kernel
    @test_throws ArgumentError Peridynamics.get_q_dim(:UNKNOWN)
    @test_throws ArgumentError Peridynamics.get_monomial_vector(:UNKNOWN, ΔX)
    @test_throws ArgumentError Peridynamics.get_gradient_extraction_matrix(:UNKNOWN)
end

@testitem "reproducing kernel mathematical properties" begin
    using Peridynamics.StaticArrays, Peridynamics.LinearAlgebra

    # Test that gradient matrices have correct properties
    monomials = [:C1, :RK1, :RK2, :PD2]

    for kernel in monomials
        Q∇ᵀ = Peridynamics.get_gradient_extraction_matrix(kernel)
        q_dim = Peridynamics.get_q_dim(kernel)

        # Check dimensions
        @test size(Q∇ᵀ) == (3, q_dim)

        # For the gradient operator, check that it correctly extracts gradients
        # For example, for RK1 with basis [1, x, y, z], the gradient should pick out [x, y, z] derivatives
        if kernel == :RK1
            # Q∇ᵀ should extract [∂/∂x, ∂/∂y, ∂/∂z] from [1, x, y, z]
            @test Q∇ᵀ[1, 2] == 1  # ∂x/∂x = 1
            @test Q∇ᵀ[2, 3] == 1  # ∂y/∂y = 1
            @test Q∇ᵀ[3, 4] == 1  # ∂z/∂z = 1
            @test Q∇ᵀ[1, 1] == 0  # ∂1/∂x = 0
        elseif kernel == :C1
            # For C1 with basis [x, y, z], gradient should be identity
            @test Q∇ᵀ == I
        elseif kernel == :PD2
            # For PD2 with basis [x, y, z, x², xy, xz, y², yz, z²]
            # Gradient should pick out linear terms
            @test Q∇ᵀ[1, 1] == 1  # ∂x/∂x = 1
            @test Q∇ᵀ[2, 2] == 1  # ∂y/∂y = 1
            @test Q∇ᵀ[3, 3] == 1  # ∂z/∂z = 1
        end
    end
end

@testitem "reproducing kernel material constructor validation" begin
    # Test that all kernels work with material constructor
    monomials = [:C1, :RK1, :RK2, :PD2]

    for kernel in monomials
        mat = RKCMaterial(monomial=kernel)
        @test mat.monomial == kernel
        @test Peridynamics.get_q_dim(kernel) > 0
    end

    # Test invalid kernel
    @test_throws ArgumentError RKCMaterial(monomial=:INVALID)
end

@testitem "reproducing kernel basic functionality test" begin
    using Peridynamics.StaticArrays, Peridynamics.LinearAlgebra

    # Test each kernel with a simple configuration
    monomials = [:C1, :RK1, :RK2, :PD2]

    for monomial in monomials
        # Create a simple test system with the kernel
        pos, vol = uniform_box(1, 1, 1, 0.4)
        body = Body(RKCMaterial(; monomial), pos, vol)
        material!(body; horizon=1.5, rho=1, E=210e9, nu=0.25, Gc=1.0)

        dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
        chunk = dh.chunks[1]
        (; storage) = chunk

        # Test that the system initializes without errors
        @test length(storage.defgrad) > 0
        @test length(storage.gradient_weight) > 0

        # Test that calc_force_density runs without errors
        @test_nowarn Peridynamics.calc_force_density!(dh, 0.0, 0.0)

        # Test that deformation gradients are computed (non-zero)
        F_total = sum(abs.(storage.defgrad))
        @test F_total > 0.0
    end
end

@testitem "Affected points RKCMaterial NewtonRaphson" begin
    using Peridynamics: NewtonRaphson
    position = [0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    volume = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    body = Body(RKCMaterial(), position, volume)
    material!(body; horizon=1.1, rho=1, E=210e9, nu=0.25)
    ts = NewtonRaphson(steps=1)
    dh = Peridynamics.threads_data_handler(body, ts, 1)
    chunk = dh.chunks[1]
    (; storage, system) = chunk
    (; affected_points) = storage

    @test Peridynamics.get_affected_points(system, 1) == [1, 2]
    @test affected_points[1] == [1, 2, 3] # one layer of neighbors more

    @test Peridynamics.get_affected_points(system, 2) == [1, 2, 3]
    @test affected_points[2] == [1, 2, 3, 4]

    @test Peridynamics.get_affected_points(system, 3) == [2, 3, 4]
    @test affected_points[3] == [1, 2, 3, 4, 5]

    @test Peridynamics.get_affected_points(system, 4) == [3, 4, 5]
    @test affected_points[4] == [2, 3, 4, 5, 6]

    @test Peridynamics.get_affected_points(system, 5) == [4, 5, 6]
    @test affected_points[5] == [3, 4, 5, 6, 7]

    @test Peridynamics.get_affected_points(system, 6) == [5, 6, 7]
    @test affected_points[6] == [4, 5, 6, 7, 8]

    @test Peridynamics.get_affected_points(system, 7) == [6, 7, 8]
    @test affected_points[7] == [5, 6, 7, 8, 9]

    @test Peridynamics.get_affected_points(system, 8) == [7, 8, 9]
    @test affected_points[8] == [6, 7, 8, 9, 10]

    @test Peridynamics.get_affected_points(system, 9) == [8, 9, 10]
    @test affected_points[9] == [7, 8, 9, 10]

    @test Peridynamics.get_affected_points(system, 10) == [9, 10]
    @test affected_points[10] == [8, 9, 10]
end
