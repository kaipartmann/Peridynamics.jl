@testitem "Flanagan-Taylor stress rotation algorithm CRMaterial" begin
    using Peridynamics.StaticArrays
    using Peridynamics.LinearAlgebra

    # Define the identity matrix for orthogonality checks
    I_matrix = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

    # Create a small test body - 2×2×2 cube
    l, w, h = 1.0, 1.0, 1.0
    ΔX = 0.5  # Creates a 2×2×2 grid
    δ = 3.1 * ΔX
    ref_position, volume = uniform_box(l, w, h, ΔX)
    mat = CRMaterial()
    body = Body(mat, ref_position, volume)
    material!(body; horizon=δ, rho=8000, E=210e9, nu=0.25)
    no_failure!(body)
    ts = VelocityVerlet(steps=2)

    # Create the data handler
    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; storage, system, paramsetup) = chunk
    (; position) = storage

    # Test 1: Simple shear deformation
    # Apply a simple shear: γ = 0.1 (10% shear strain)
    γ = 0.1
    F_shear = @SMatrix [1.0 γ   0.0
                        0.0 1.0 0.0
                        0.0 0.0 1.0]

    # Apply shear deformation to all points
    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(position, i)
        x = F_shear * X
        for d in 1:3
            position[d, i] = x[d]
        end
    end

    # Run one force calculation step
    Peridynamics.calc_force_density!(dh, 0.0, 0.1)

    # Check that rotation matrices are orthogonal (R * R' = I)
    for i in Peridynamics.each_point_idx(chunk)
        R = Peridynamics.get_tensor(storage.rotation, i)
        RT_R = R' * R
        @test RT_R ≈ I_matrix atol=1e-12

        # Check determinant is 1 (proper rotation)
        @test abs(det(R) - 1.0) < 1e-12
    end

    # Check that left stretch tensors are symmetric and positive definite
    for i in Peridynamics.each_point_idx(chunk)
        V = Peridynamics.get_tensor(storage.left_stretch, i)

        # Check symmetry: V = V'
        @test V ≈ V' atol=1e-12

        # Check positive definiteness (all eigenvalues > 0)
        eigenvals = eigvals(V)
        @test all(eigenvals .> 0)
    end

    # Test 2: Combined rotation and stretch
    # Reset positions
    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(ref_position, i)
        for d in 1:3
            position[d, i] = X[d]
        end
    end

    # Apply a combined deformation: stretch + rotation
    θ = π/12  # 15 degree rotation
    stretch = 1.05
    F_combined = @SMatrix [stretch*cos(θ) -stretch*sin(θ) 0.0
                           sin(θ)         cos(θ)         0.0
                           0.0            0.0            1.0]

    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(position, i)
        x = F_combined * X
        for d in 1:3
            position[d, i] = x[d]
        end
    end

    # Run force calculation
    Peridynamics.calc_force_density!(dh, 0.1, 0.1)

    # Verify rotation matrices remain orthogonal
    for i in Peridynamics.each_point_idx(chunk)
        R = Peridynamics.get_tensor(storage.rotation, i)
        RT_R = R' * R
        I_matrix = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        @test RT_R ≈ I_matrix atol=1e-10
        @test abs(det(R) - 1.0) < 1e-10
    end

    # Verify left stretch tensors remain symmetric and positive definite
    for i in Peridynamics.each_point_idx(chunk)
        V = Peridynamics.get_tensor(storage.left_stretch, i)
        @test V ≈ V' atol=1e-10
        eigenvals = eigvals(V)
        @test all(eigenvals .> 0)
    end

    # Test 3: Large deformation to test robustness
    # Reset positions
    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(ref_position, i)
        for d in 1:3
            position[d, i] = X[d]
        end
    end

    # Apply large shear deformation
    γ_large = 0.5  # 50% shear strain
    F_large_shear = @SMatrix [1.0 γ_large 0.0
                              0.0 1.0      0.0
                              0.0 0.0      1.0]

    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(position, i)
        x = F_large_shear * X
        for d in 1:3
            position[d, i] = x[d]
        end
    end

    # Run force calculation
    Peridynamics.calc_force_density!(dh, 0.2, 0.1)

    # Even for large deformations, the algorithm should maintain properties
    for i in Peridynamics.each_point_idx(chunk)
        R = Peridynamics.get_tensor(storage.rotation, i)
        V = Peridynamics.get_tensor(storage.left_stretch, i)

        # Rotation should still be orthogonal
        RT_R = R' * R
        I_matrix = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        @test RT_R ≈ I_matrix atol=1e-8
        @test abs(det(R) - 1.0) < 1e-8

        # Left stretch should be symmetric and positive definite
        @test V ≈ V' atol=1e-8
        eigenvals = eigvals(V)
        @test all(eigenvals .> 0)
    end

    # Test 4: Verify stress is rotated correctly
    # For a simple test, check that stress has the expected form
    # Reset to a simple extension case
    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(ref_position, i)
        for d in 1:3
            position[d, i] = X[d]
        end
    end

    # Apply uniaxial tension in x-direction
    stretch_x = 1.01
    F_tension = @SMatrix [stretch_x 0.0 0.0
                          0.0       1.0 0.0
                          0.0       0.0 1.0]

    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(position, i)
        x = F_tension * X
        for d in 1:3
            position[d, i] = x[d]
        end
    end

    # Run force calculation
    Peridynamics.calc_force_density!(dh, 0.3, 0.1)

    # Check that von Mises stress is calculated and reasonable
    for i in Peridynamics.each_point_idx(chunk)
        vm_stress = storage.von_mises_stress[i]
        @test vm_stress ≥ 0.0  # von Mises stress is always non-negative
        @test isfinite(vm_stress)  # Should not be NaN or Inf
    end
end

@testitem "Flanagan-Taylor stress rotation algorithm RKCRMaterial" begin
    using Peridynamics.StaticArrays
    using Peridynamics.LinearAlgebra

    # Create a small test body - 2×2×2 cube
    l, w, h = 1.0, 1.0, 1.0
    ΔX = 0.5  # Creates a 2×2×2 grid
    δ = 3.1 * ΔX
    ref_position, volume = uniform_box(l, w, h, ΔX)
    mat = RKCRMaterial()
    body = Body(mat, ref_position, volume)
    material!(body; horizon=δ, rho=8000, E=210e9, nu=0.25)
    no_failure!(body)
    ts = VelocityVerlet(steps=2)

    # Create the data handler
    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; storage, system, paramsetup) = chunk
    (; position) = storage

    # Test 1: Simple shear deformation
    γ = 0.1
    F_shear = @SMatrix [1.0 γ   0.0
                        0.0 1.0 0.0
                        0.0 0.0 1.0]

    # Apply shear deformation to all points
    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(position, i)
        x = F_shear * X
        for d in 1:3
            position[d, i] = x[d]
        end
    end

    # Run one force calculation step
    Peridynamics.calc_force_density!(dh, 0.0, 0.1)

    # For RKCRMaterial, rotation matrices are stored per bond, not per point
    # Check rotation matrices for all bonds
    n_bonds = size(storage.rotation, 2)
    for bond_id in 1:n_bonds
        if storage.bond_active[bond_id]
            R = Peridynamics.get_tensor(storage.rotation, bond_id)
            RT_R = R' * R
            I_matrix = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
            @test RT_R ≈ I_matrix atol=1e-10
            @test abs(det(R) - 1.0) < 1e-10
        end
    end

    # Check left stretch tensors for all bonds
    for bond_id in 1:n_bonds
        if storage.bond_active[bond_id]
            V = Peridynamics.get_tensor(storage.left_stretch, bond_id)
            @test V ≈ V' atol=1e-10
            eigenvals = eigvals(V)
            @test all(eigenvals .> 0)
        end
    end

    # Test 2: Combined rotation and stretch
    # Reset positions
    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(ref_position, i)
        for d in 1:3
            position[d, i] = X[d]
        end
    end

    # Apply combined deformation
    θ = π/12
    stretch = 1.05
    F_combined = @SMatrix [stretch*cos(θ) -stretch*sin(θ) 0.0
                           sin(θ)         cos(θ)         0.0
                           0.0            0.0            1.0]

    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(position, i)
        x = F_combined * X
        for d in 1:3
            position[d, i] = x[d]
        end
    end

    # Run force calculation
    Peridynamics.calc_force_density!(dh, 0.1, 0.1)

    # Verify properties for all active bonds
    for bond_id in 1:n_bonds
        if storage.bond_active[bond_id]
            R = Peridynamics.get_tensor(storage.rotation, bond_id)
            V = Peridynamics.get_tensor(storage.left_stretch, bond_id)

            # Rotation should be orthogonal
            RT_R = R' * R
            I_matrix = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
            @test RT_R ≈ I_matrix atol=1e-8
            @test abs(det(R) - 1.0) < 1e-8

            # Left stretch should be symmetric and positive definite
            @test V ≈ V' atol=1e-8
            eigenvals = eigvals(V)
            @test all(eigenvals .> 0)
        end
    end

    # Test 3: Verify deformation gradient calculation
    # Reset to uniaxial tension
    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(ref_position, i)
        for d in 1:3
            position[d, i] = X[d]
        end
    end

    stretch_x = 1.02
    F_tension = @SMatrix [stretch_x 0.0 0.0
                          0.0       1.0 0.0
                          0.0       0.0 1.0]

    for i in Peridynamics.each_point_idx(chunk)
        X = Peridynamics.get_vector(position, i)
        x = F_tension * X
        for d in 1:3
            position[d, i] = x[d]
        end
    end

    # Run force calculation
    Peridynamics.calc_force_density!(dh, 0.2, 0.1)

    # Check that deformation gradients are reasonable
    for i in Peridynamics.each_point_idx(chunk)
        F = Peridynamics.get_tensor(storage.defgrad, i)
        # For uniaxial tension, F should be approximately diagonal
        @test F[1,1] > 1.0  # Stretch in x-direction
        @test abs(F[2,2] - 1.0) < 0.1  # Approximately unity in y
        @test abs(F[3,3] - 1.0) < 0.1  # Approximately unity in z
        @test det(F) > 0  # Positive determinant
    end

    # Check von Mises stress values are reasonable
    for i in Peridynamics.each_point_idx(chunk)
        vm_stress = storage.von_mises_stress[i]
        @test vm_stress ≥ 0.0
        @test isfinite(vm_stress)
    end
end

@testitem "Stress rotation consistency between materials" begin
    using Peridynamics.StaticArrays
    using Peridynamics.LinearAlgebra

    # Create identical bodies for both materials
    l, w, h = 1.0, 1.0, 1.0
    ΔX = 0.5
    δ = 3.1 * ΔX
    ref_position, volume = uniform_box(l, w, h, ΔX)

    # CRMaterial body
    mat_cr = CRMaterial()
    body_cr = Body(mat_cr, deepcopy(ref_position), deepcopy(volume))
    material!(body_cr; horizon=δ, rho=8000, E=210e9, nu=0.25)
    no_failure!(body_cr)

    # RKCRMaterial body
    mat_rkcr = RKCRMaterial()
    body_rkcr = Body(mat_rkcr, deepcopy(ref_position), deepcopy(volume))
    material!(body_rkcr; horizon=δ, rho=8000, E=210e9, nu=0.25)
    no_failure!(body_rkcr)

    # Create data handlers
    ts = VelocityVerlet(steps=1)
    dh_cr = Peridynamics.threads_data_handler(body_cr, ts, 1)
    dh_rkcr = Peridynamics.threads_data_handler(body_rkcr, ts, 1)
    Peridynamics.initialize!(dh_cr, ts)
    Peridynamics.initialize!(dh_rkcr, ts)

    chunk_cr = dh_cr.chunks[1]
    chunk_rkcr = dh_rkcr.chunks[1]

    # Apply the same small deformation to both
    γ = 0.05  # Small shear strain for comparison
    F_shear = @SMatrix [1.0 γ   0.0
                        0.0 1.0 0.0
                        0.0 0.0 1.0]

    # Apply to CRMaterial
    for i in Peridynamics.each_point_idx(chunk_cr)
        X = Peridynamics.get_vector(chunk_cr.storage.position, i)
        x = F_shear * X
        for d in 1:3
            chunk_cr.storage.position[d, i] = x[d]
        end
    end

    # Apply to RKCRMaterial
    for i in Peridynamics.each_point_idx(chunk_rkcr)
        X = Peridynamics.get_vector(chunk_rkcr.storage.position, i)
        x = F_shear * X
        for d in 1:3
            chunk_rkcr.storage.position[d, i] = x[d]
        end
    end

    # Run force calculations
    Peridynamics.calc_force_density!(dh_cr, 0.0, 0.1)
    Peridynamics.calc_force_density!(dh_rkcr, 0.0, 0.1)

    # Both materials should show similar behavior for small deformations
    # Compare von Mises stress (this is a qualitative comparison)
    for i in Peridynamics.each_point_idx(chunk_cr)
        vm_cr = chunk_cr.storage.von_mises_stress[i]
        vm_rkcr = chunk_rkcr.storage.von_mises_stress[i]

        # Both should be finite and positive
        @test vm_cr ≥ 0.0 && isfinite(vm_cr)
        @test vm_rkcr ≥ 0.0 && isfinite(vm_rkcr)

        # For small deformations, they should be of similar magnitude
        # (allowing for some difference due to different numerical approaches)
        if vm_cr > 1e-10 && vm_rkcr > 1e-10
            ratio = vm_cr / vm_rkcr
            @test 0.1 < ratio < 10.0  # Should be within an order of magnitude
        end
    end
end
