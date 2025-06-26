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
    @test mat1.kernel == cubic_b_spline_kernel
    @test mat1.constitutive_model isa LinearElastic
    @test mat1.dmgmodel isa CriticalStretch
    @test mat1.maxdmg == 0.85
    @test mat1.reprkernel == :C1
    @test mat1.regfactor == 1e-13

    # Test constructor with parameters
    mat2 = RKCMaterial(
        kernel = linear_kernel,
        model = SaintVenantKirchhoff(),
        dmgmodel = CriticalStretch(),
        maxdmg = 0.75,
        reprkernel = :C1,
        regfactor = 1e-10
    )
    @test mat2.kernel == linear_kernel
    @test mat2.constitutive_model isa SaintVenantKirchhoff
    @test mat2.dmgmodel isa CriticalStretch
    @test mat2.maxdmg == 0.75
    @test mat2.reprkernel == :C1
    @test mat2.regfactor == 1e-10

    # Test constructor with invalid regfactor
    @test_throws ArgumentError RKCMaterial(regfactor = 2.0)
    @test_throws ArgumentError RKCMaterial(regfactor = -0.5)
end

@testitem "RKCRMaterial initialization" begin
    # Test default constructor
    mat1 = RKCRMaterial()
    @test mat1.kernel == cubic_b_spline_kernel
    @test mat1.constitutive_model isa LinearElastic
    @test mat1.dmgmodel isa CriticalStretch
    @test mat1.maxdmg == 0.85
    @test mat1.reprkernel == :C1
    @test mat1.regfactor == 1e-13

    # Test constructor with parameters (only LinearElastic is supported)
    mat2 = RKCRMaterial(
        kernel = linear_kernel,
        model = LinearElastic(),
        dmgmodel = CriticalStretch(),
        maxdmg = 0.75,
        reprkernel = :C1,
        regfactor = 1e-10
    )
    @test mat2.kernel == linear_kernel
    @test mat2.constitutive_model isa LinearElastic
    @test mat2.dmgmodel isa CriticalStretch
    @test mat2.maxdmg == 0.75
    @test mat2.reprkernel == :C1
    @test mat2.regfactor == 1e-10

    # Test failure with non-LinearElastic model
    @test_throws ArgumentError RKCRMaterial(model = SaintVenantKirchhoff())

    # Test constructor with invalid regfactor
    @test_throws ArgumentError RKCRMaterial(regfactor = 2.0)
    @test_throws ArgumentError RKCRMaterial(regfactor = -0.5)
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

# @testitem "consistency between RKCMaterial and RKCRMaterial for small deformation" begin
#     pos, vol = uniform_box(1, 1, 1, 0.1)

#     body1 = Body(RKCMaterial(), pos, vol)
#     material!(body1; horizon=1.5, rho=1, E=210e9, nu=0.25, Gc=1.0)

#     body2 = Body(RKCRMaterial(), deepcopy(pos), deepcopy(vol))
#     material!(body2; horizon=1.5, rho=1, E=210e9, nu=0.25, Gc=1.0)

#     dh1 = Peridynamics.threads_data_handler(body1, VelocityVerlet(steps=1), 1)
#     chunk1 = dh1.chunks[1]

#     dh2 = Peridynamics.threads_data_handler(body2, VelocityVerlet(steps=1), 1)
#     chunk2 = dh2.chunks[1]

#     # Apply small deformation
#     stretch_factor = 1.00001
#     for i in eachindex(vol)
#         chunk1.storage.displacement[:,i] = [chunk1.storage.position[1,i] * (stretch_factor - 1), 0.0, 0.0]
#         chunk1.storage.position[:,i] += chunk1.storage.displacement[:,i]

#         chunk2.storage.displacement[:,i] = [chunk2.storage.position[1,i] * (stretch_factor - 1), 0.0, 0.0]
#         chunk2.storage.position[:,i] += chunk2.storage.displacement[:,i]
#     end

#     # Calculate force density
#     Peridynamics.calc_force_density!(dh1, 0.0, 0.0)
#     Peridynamics.calc_force_density!(dh2, 0.0, 0.0)

#     # For small deformations, internal forces should be very similar
#     @test chunk1.storage.b_int ≈ chunk2.storage.b_int rtol=1e-4
# end
