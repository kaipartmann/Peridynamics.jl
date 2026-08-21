@testitem "Piola transform of cauchy stress" begin
    # deformation gradient
    F = [2.0 1.0 0.0; 1.0 2.0 1.0; 1.0 1.0 2.0]

    # calc σ from P and back to P
    P = [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0]
    σ = Peridynamics.cauchy_stress(P, F)
    @test Peridynamics.first_piola_kirchhoff(σ, F) ≈ P

    # calc P from σ and back to σ
    σ = [1.0 2.0 3.0; 2.0 4.0 5.0; 3.0 5.0 6.0]
    P = Peridynamics.first_piola_kirchhoff(σ, F)
    @test Peridynamics.cauchy_stress(P, F) ≈ σ
end

@testitem "Von Mises stress calculation" begin
    using Peridynamics.StaticArrays

    # zero stress
    σx, σy, σz = 0, 0, 0
    τxy, τyz, τzx = 0, 0, 0
    σ = @SMatrix [σx τxy τzx; τxy σy τyz; τzx τyz σz]
    σvm = Peridynamics.von_mises_stress(σ)
    @test σvm ≈ 0 atol=eps()

    # uniaxial stress x
    σx, σy, σz = 100.0, 0, 0
    τxy, τyz, τzx = 0, 0, 0
    σ = @SMatrix [σx τxy τzx; τxy σy τyz; τzx τyz σz]
    σvm = Peridynamics.von_mises_stress(σ)
    @test σvm ≈ σx

    # uniaxial stress y
    σx, σy, σz = 0, 100.0, 0
    τxy, τyz, τzx = 0, 0, 0
    σ = @SMatrix [σx τxy τzx; τxy σy τyz; τzx τyz σz]
    σvm = Peridynamics.von_mises_stress(σ)
    @test σvm ≈ σy

    # uniaxial stress z
    σx, σy, σz = 0, 0, 100.0
    τxy, τyz, τzx = 0, 0, 0
    σ = @SMatrix [σx τxy τzx; τxy σy τyz; τzx τyz σz]
    σvm = Peridynamics.von_mises_stress(σ)
    @test σvm ≈ σz

    # shear and axial loading
    σx, σy, σz = 500, 300, 400
    τxy, τyz, τzx = 300, 0, 0
    σ = @SMatrix [σx τxy τzx; τxy σy τyz; τzx τyz σz]
    σvm = Peridynamics.von_mises_stress(σ)
    @test σvm ≈ √300000

    # shear and axial loading
    σx, σy, σz = 500, 300, 400
    τxy, τyz, τzx = 0, 300, 0
    σ = @SMatrix [σx τxy τzx; τxy σy τyz; τzx τyz σz]
    σvm = Peridynamics.von_mises_stress(σ)
    @test σvm ≈ √300000

    # shear and axial loading
    σx, σy, σz = 500, 300, 400
    τxy, τyz, τzx = 0, 0, 300
    σ = @SMatrix [σx τxy τzx; τxy σy τyz; τzx τyz σz]
    σvm = Peridynamics.von_mises_stress(σ)
    @test σvm ≈ √300000

    # all shear and axial loading
    σx, σy, σz = 100, 100, 100
    τxy, τyz, τzx = 100, 100, 100
    σ = @SMatrix [σx τxy τzx; τxy σy τyz; τzx τyz σz]
    σvm = Peridynamics.von_mises_stress(σ)
    @test σvm ≈ 300
end

@testitem "init_stress_rotation!" begin
    using Peridynamics.StaticArrays

    ref_position = [0.0 1.0 0.0 0.0 2.0
                    0.0 0.0 1.0 0.0 2.0
                    0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    δ = 1.5
    body = Body(CRMaterial(model=LinearElastic()), ref_position, volume)
    material!(body, horizon=δ, rho=1, E=1, nu=0.25)
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    (; mat, storage, system, paramsetup) = dh.chunks[1]

    F = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    Ḟ = @SMatrix [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    Δt = 0.1
    i = 1

    Dᵣ = Peridynamics.init_stress_rotation!(storage, F, Ḟ, Δt, i)
    D_solution = @SMatrix [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    @test Dᵣ ≈ D_solution atol=10eps()

    F = @SMatrix [1.1 0.1 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    Ḟ = @SMatrix [0.1 0.0 0.0; 0.0 0.1 0.0; 0.0 0.0 0.1]
    Δt = 0.1

    Dᵣ = Peridynamics.init_stress_rotation!(storage, F, Ḟ, Δt, i)
    D_solution = @SMatrix [0.09090499812142767 -0.004541322313765097 0.0;
                           -0.004541322313765098 0.10000417167571606 0.0;
                           0.0 0.0 0.1]
    @test Dᵣ ≈ D_solution
end

# ---- the Flanagan-Taylor stress rotation used by the rotated formulations ----

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
        @test RT_R ≈ I atol=1e-10
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
        @test RT_R ≈ I atol=1e-8
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

@testitem "init_stress_rotation!: a singular deformation gradient resets the rotation" setup=[Fixtures] begin
    using Peridynamics.StaticArrays
    body = Fixtures.tetra4(CRMaterial(model=LinearElastic()))
    c = Fixtures.chunk(body)
    (; storage) = c
    F = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    Ḟ = @SMatrix [0.1 0.0 0.0; 0.0 0.1 0.0; 0.0 0.0 0.1]
    D = Peridynamics.init_stress_rotation!(storage, F, Ḟ, 0.1, 1)
    @test D ≈ Ḟ
    @test Peridynamics.get_tensor(storage.rotation, 1) ≈ F
    # F⁻¹ contains NaN for a singular F: the rotation and the stretch are zeroed, D is zero
    F = @SMatrix [1.0 0.0 0.0; 1.0 0.0 0.0; 1.0 0.0 0.0]
    D = Peridynamics.init_stress_rotation!(storage, F, Ḟ, 0.1, 1)
    @test iszero(D)
    @test iszero(Peridynamics.get_tensor(storage.rotation, 1))
    @test iszero(Peridynamics.get_tensor(storage.left_stretch, 1))
end

@testitem "sym_eigvals: closed-form eigenvalues of a symmetric tensor" begin
    using Peridynamics.StaticArrays
    using Peridynamics.LinearAlgebra

    sym_eigvals = Peridynamics.sym_eigvals

    # a diagonal tensor: the spectrum is the diagonal, and it comes back sorted descending
    @test sym_eigvals(SMatrix{3,3,Float64,9}(Diagonal([1.0, 9.0, 4.0]))) ==
          SVector{3,Float64}(9.0, 4.0, 1.0)
    @test sym_eigvals(SMatrix{3,3,Float64,9}(Diagonal([1.0, 4.0, 1.0]))) ==
          SVector{3,Float64}(4.0, 1.0, 1.0)

    # repeated eigenvalues
    @test sym_eigvals(2.5 * SMatrix{3,3,Float64,9}(I)) ≈ SVector{3,Float64}(2.5, 2.5, 2.5)

    # a general symmetric tensor, against a full eigendecomposition
    A = SMatrix{3,3,Float64,9}(4.0, 1.0, 0.7, 1.0, 3.0, -0.4, 0.7, -0.4, 2.0)
    @test sym_eigvals(A) ≈ SVector{3,Float64}(reverse(eigvals(Symmetric(Matrix(A)))))

    # the invariants are reproduced exactly enough
    λ = sym_eigvals(A)
    @test sum(λ) ≈ tr(A)
    @test prod(λ) ≈ det(A)

    # always sorted descending
    @test issorted(sym_eigvals(A); rev=true)
end

@testitem "hencky_and_invstretch: logarithmic strain and inverse right stretch" begin
    using Peridynamics.StaticArrays
    using Peridynamics.LinearAlgebra

    hencky = Peridynamics.hencky_and_invstretch
    Id = SMatrix{3,3,Float64,9}(I)

    # the reference, built from a full eigendecomposition of C
    function reference(C)
        E = eigen(Symmetric(Matrix(C)))
        ε = zeros(3, 3)
        Uinv = zeros(3, 3)
        for k in 1:3
            n = E.vectors[:, k]
            ε .+= 0.5 * log(E.values[k]) .* (n * n')
            Uinv .+= 1 / sqrt(E.values[k]) .* (n * n')
        end
        return SMatrix{3,3}(ε), SMatrix{3,3}(Uinv)
    end

    # no deformation
    let (ε, Uinv) = hencky(Id)
        @test ε ≈ zero(Id) atol=1e-14
        @test Uinv ≈ Id
    end

    testcases = (
        # a triple eigenvalue
        1.44 * Id,
        # uniaxial stretch: a double eigenvalue on the small pair
        SMatrix{3,3,Float64,9}(Diagonal([1.4, 1.0, 1.0]))' *
        SMatrix{3,3,Float64,9}(Diagonal([1.4, 1.0, 1.0])),
        # equibiaxial: a double eigenvalue on the large pair
        SMatrix{3,3,Float64,9}(Diagonal([1.3, 1.3, 1.0]))' *
        SMatrix{3,3,Float64,9}(Diagonal([1.3, 1.3, 1.0])),
        # an unsorted diagonal, which must not be mistaken for three distinct eigenvalues
        SMatrix{3,3,Float64,9}(Diagonal([1.0, 4.0, 1.0])),
        # three distinct eigenvalues, with shear
        let F = SMatrix{3,3,Float64,9}(1.2, 0.1, 0.05, -0.07, 0.95, 0.02, 0.03, -0.04, 1.05)
            F' * F
        end,
        # two eigenvalues that are close but not equal
        SMatrix{3,3,Float64,9}(Diagonal([2.0, 1.0 + 1e-10, 1.0])),
    )

    for C in testcases
        εc, Uinvc = hencky(C)
        εref, Uref = reference(C)
        @test εc ≈ εref atol=1e-9
        @test Uinvc ≈ Uref atol=1e-9
        # U⁻¹ is the inverse square root of C, whatever the multiplicity of its eigenvalues
        @test Uinvc * C * Uinvc ≈ Id atol=1e-9
        # both are symmetric, because they are isotropic functions of a symmetric tensor
        @test εc ≈ εc' atol=1e-12
        @test Uinvc ≈ Uinvc' atol=1e-12
    end

    # uniaxial stretch has the closed form ε₁₁ = ln λ
    let λ = 1.4, (εu, _) = hencky(SMatrix{3,3,Float64,9}(Diagonal([λ^2, 1.0, 1.0])))
        @test εu[1, 1] ≈ log(λ)
        @test εu[2, 2] ≈ 0.0 atol=1e-14
        @test εu[3, 3] ≈ 0.0 atol=1e-14
    end

    # an isochoric deformation has a traceless Hencky strain
    let F = SMatrix{3,3,Float64,9}(Diagonal([2.0, 1 / sqrt(2), 1 / sqrt(2)]))
        @test tr(first(hencky(F' * F))) ≈ 0.0 atol=1e-14
    end

    # the float type of the input is carried through
    C32 = SMatrix{3,3,Float32,9}(Diagonal([1.44f0, 1.0f0, 0.81f0]))
    ε32, U32 = hencky(C32)
    @test eltype(ε32) === Float32
    @test eltype(U32) === Float32

    # nothing is allocated, which is what the closed form is for
    C = SMatrix{3,3,Float64,9}(1.44, 0.1, 0.0, 0.1, 1.0, 0.0, 0.0, 0.0, 0.81)
    hencky(C)
    @test (@allocated hencky(C)) == 0
end
