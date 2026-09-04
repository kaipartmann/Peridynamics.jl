# `RKCMaterial`: the reproducing kernel correspondence formulation of
# `src/physics/rk_correspondence.jl`.

@testitem "damage changed flag" begin
    pos, vol = uniform_box(1, 1, 1, 0.4)
    body = Body(RKCMaterial(kernel=cubic_b_spline_kernel_norm), pos, vol)
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
    @test mat1.kernel == const_one_kernel
    @test mat1.constitutive_model isa SaintVenantKirchhoff
    @test mat1.dmgmodel isa CriticalStretch
    @test Peridynamics.monomial(mat1) == :C1
    @test mat1.epsilon == 1e-3
    @test mat1.lambda == 0
    @test mat1.beta == 0

    # Test constructor with parameters, `lambda` or `beta` select the legacy regularization
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
    @test Peridynamics.monomial(mat2) == :C1
    @test mat2.epsilon == 0
    @test mat2.lambda == 0
    @test mat2.beta == 1e-8

    # Setting only one of the legacy parameters keeps the legacy default of the other one
    mat3 = RKCMaterial(lambda = 1e-6)
    @test mat3.epsilon == 0
    @test mat3.lambda == 1e-6
    @test mat3.beta ≈ sqrt(eps())

    # Test constructor with a singular value floor
    mat4 = RKCMaterial(epsilon = 1e-2)
    @test mat4.epsilon == 1e-2
    @test mat4.lambda == 0
    @test mat4.beta == 0

    # `epsilon = 0` is the exact pseudo-inverse, not the legacy regularization
    mat5 = RKCMaterial(epsilon = 0)
    @test mat5.epsilon == 0
    @test mat5.lambda == 0
    @test mat5.beta == 0

    # Test constructor with invalid epsilon/lambda/beta
    @test_throws ArgumentError RKCMaterial(epsilon = -0.5)
    @test_throws ArgumentError RKCMaterial(lambda = -0.5)
    @test_throws ArgumentError RKCMaterial(beta = -0.5)

    # The two regularizations cannot be combined
    @test_throws ArgumentError RKCMaterial(epsilon = 1e-3, lambda = 1e-6)
    @test_throws ArgumentError RKCMaterial(epsilon = 1e-3, beta = 1e-6)
end

@testitem "RKCMaterial: the monomial is a type parameter" begin
    # the monomial basis selects the size of the moment matrix, so it is part of the type
    # and `monomial(mat)` is a compile-time constant
    for monomial in (:C1, :RK1, :RK2, :PD2)
        mat = RKCMaterial(; monomial)
        @test mat isa RKCMaterial{SaintVenantKirchhoff,typeof(const_one_kernel),
                                  CriticalStretch,monomial}
        @test mat isa Peridynamics.AbstractRKCMaterial{SaintVenantKirchhoff,
                                                       Peridynamics.NoCorrection,monomial}
        @test Peridynamics.monomial(mat) == monomial
        @test (@inferred Peridynamics.monomial(mat)) == monomial
        @test !hasproperty(mat, :monomial)
        matr = RKCRMaterial(; monomial)
        @test matr isa RKCRMaterial{SaintVenantKirchhoff,typeof(const_one_kernel),
                                    CriticalStretch,monomial}
        @test Peridynamics.monomial(matr) == monomial
    end
    @test RKCMaterial(monomial=:C1) !== RKCMaterial(monomial=:RK1)
    @test typeof(RKCMaterial(monomial=:C1)) != typeof(RKCMaterial(monomial=:RK1))
end

@testitem "get_invreg_params: exactly one regularization is active" begin
    import Peridynamics: get_invreg_params

    # nothing given: the adaptive regularization with its default floor
    @test get_invreg_params(nothing, nothing, nothing) == (1e-3, 0.0, 0.0)
    @test get_invreg_params(1e-2, nothing, nothing) == (1e-2, 0.0, 0.0)
    @test get_invreg_params(0, nothing, nothing) == (0.0, 0.0, 0.0)
    @test get_invreg_params(1, nothing, nothing) === (1.0, 0.0, 0.0) # converted to Float64

    # a legacy keyword switches the floor off and fills in the legacy default of the other
    @test get_invreg_params(nothing, 1e-6, nothing) == (0.0, 1e-6, sqrt(eps()))
    @test get_invreg_params(nothing, nothing, 1e-8) == (0.0, 0.0, 1e-8)
    @test get_invreg_params(nothing, 1e-6, 1e-8) == (0.0, 1e-6, 1e-8)
    @test get_invreg_params(nothing, 0, 0) === (0.0, 0.0, 0.0)

    # negative values and mixing the two regularizations are rejected
    @test_throws ArgumentError get_invreg_params(-1e-3, nothing, nothing)
    @test_throws ArgumentError get_invreg_params(nothing, -1e-6, nothing)
    @test_throws ArgumentError get_invreg_params(nothing, nothing, -1e-8)
    @test_throws ArgumentError get_invreg_params(1e-3, 1e-6, nothing)
    @test_throws ArgumentError get_invreg_params(1e-3, nothing, 1e-8)
    @test_throws ArgumentError get_invreg_params(0, 0, 0)
end

@testitem "log_material: the monomial and the active regularization are logged" begin
    # the monomial is a type parameter and not a field, so it needs its own line
    msg = Peridynamics.log_material(RKCMaterial(monomial=:RK2))
    @test contains(msg, "RKCMaterial")
    @test contains(msg, "monomial type")
    @test contains(msg, "RK2")
    @test contains(msg, "singular value floor")
    @test !contains(msg, "Tikhonov")
    @test !contains(msg, "SVD truncation")

    # the legacy regularization is logged with both of its parameters and no floor
    msg = Peridynamics.log_material(RKCMaterial(lambda=1e-6))
    @test contains(msg, "Tikhonov regularization parameter")
    @test contains(msg, "SVD truncation parameter")
    @test !contains(msg, "singular value floor")

    # the rotated material logs the same way
    msg = Peridynamics.log_material(RKCRMaterial(monomial=:RK1, epsilon=1e-2))
    @test contains(msg, "RKCRMaterial")
    @test contains(msg, "RK1")
    @test contains(msg, "singular value floor")
    @test !contains(msg, "Tikhonov")
end

@testitem "gradient weights calculation" begin
    using Peridynamics.LinearAlgebra

    pos, vol = uniform_box(1, 1, 1, 0.4)
    body = Body(RKCMaterial(kernel=cubic_b_spline_kernel_norm), pos, vol)
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

@testitem "RKCMaterial: deformation gradient of a homogeneous deformation" begin
    # the reproducing kernel gradient reproduces the identity and a homogeneous stretch
    # exactly, for every point of a small cube in which every point sees every other point
    using Peridynamics.StaticArrays, Peridynamics.LinearAlgebra
    pos, vol = uniform_box(1, 1, 1, 0.25)
    body = Body(RKCMaterial(kernel=cubic_b_spline_kernel_norm), pos, vol)
    material!(body; horizon=0.76, rho=1, E=210e9, nu=0.25, Gc=1.0)
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    (; storage, system) = dh.chunks[1]
    (; position, defgrad) = storage
    # no displacement: the identity
    Peridynamics.calc_force_density!(dh, 0.0, 0.0)
    @test all(isapprox(Peridynamics.get_tensor(defgrad, i, Peridynamics.dims(system)), I; atol=1e-12) for i in eachindex(vol))
    # a small uniform stretch in x
    F_a = @SMatrix [1.00001 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    for i in eachindex(vol)
        position[:, i] = F_a * position[:, i]
    end
    Peridynamics.calc_force_density!(dh, 0.0, 0.0)
    @test all(isapprox(Peridynamics.get_tensor(defgrad, i, Peridynamics.dims(system)), F_a; atol=1e-5) for i in eachindex(vol))
end

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
        @test Peridynamics.monomial(mat) == kernel
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


@testitem "RKCMaterial: show" begin
    @test contains(sprint(show, RKCMaterial()), "RKCMaterial")
    # the monomial is visible as the last type parameter
    @test contains(sprint(show, RKCMaterial(monomial=:RK1)), "RKCMaterial")
    @test contains(sprint(show, RKCMaterial(monomial=:RK1)), ":RK1")
    @test contains(sprint(show, MIME("text/plain"), RKCMaterial()), "RKCMaterial")
end

@testitem "RKCMaterial: invalid keyword values are rejected" begin
    @test_throws ArgumentError RKCMaterial(; monomial=:NoSuchMonomial)
    @test_throws ArgumentError RKCMaterial(; epsilon=-1.0)
    @test_throws ArgumentError RKCMaterial(; lambda=-1.0)
    @test_throws ArgumentError RKCMaterial(; beta=-1.0)
    @test_throws ArgumentError RKCMaterial(; epsilon=1e-3, lambda=0.0)
    @test_throws ArgumentError RKCRMaterial(; monomial=:NoSuchMonomial)
    @test_throws ArgumentError RKCRMaterial(; epsilon=-1.0)
    @test_throws ArgumentError RKCRMaterial(; lambda=-1.0)
    @test_throws ArgumentError RKCRMaterial(; beta=-1.0)
    @test_throws ArgumentError RKCRMaterial(; epsilon=1e-3, beta=0.0)
end

@testitem "rkc_weights!: the singular value floor is inactive for intact families" setup=[Fixtures] begin
    # The default floor only damps singular values below `epsilon * σ_max`. The moment
    # matrix of an intact family is far better conditioned than that, so the gradient
    # weights are bit-identical to those of the exact inverse and to the legacy
    # regularization with no truncation: the new default changes nothing in the bulk.
    function weights(mat)
        body = Fixtures.cube(mat; n=5)
        dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
        return copy(dh.chunks[1].storage.gradient_weight)
    end
    for monomial in (:C1, :RK1, :RK2, :PD2)
        @testset "$monomial" begin
            floor = weights(RKCMaterial(; monomial))
            exact = weights(RKCMaterial(; monomial, epsilon=0))
            legacy = weights(RKCMaterial(; monomial, lambda=0, beta=0))
            @test floor == exact
            @test floor == legacy
            @test all(isfinite, floor)
        end
    end
end

@testitem "rkc_weights!: the singular value floor bounds a degenerate family" begin
    # Two points and the neighbors of the first one all on the x-axis: for `:C1` the moment
    # matrix has rank one. The floor keeps the gradient weights of the degenerate family
    # finite and bounded, and the resolved direction is still reproduced exactly.
    using Peridynamics.LinearAlgebra
    pos = [0.0 1.0 -1.0 2.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
    vol = fill(1.0, 4)
    body = Body(RKCMaterial(), pos, vol)
    # a critical stretch far above the stretch applied below, so no bond fails
    material!(body; horizon=1.5, rho=1, E=210e9, nu=0.25, epsilon_c=0.1)
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    (; storage, system) = dh.chunks[1]
    Φ = [Peridynamics.get_vector(storage.gradient_weight, b, Peridynamics.dims(system))
         for b in Peridynamics.each_bond_idx(system, 1)]
    @test all(v -> all(isfinite, v), Φ)
    # a bond along x contributes only to ∂/∂x, the singular directions get no weight
    @test all(v -> v[2] == 0 && v[3] == 0, Φ)
    # the reproducing condition Σ Φ ⊗ ΔX = I holds in the resolved direction
    ΔX = [Peridynamics.get_vector_diff(system.position, 1, Peridynamics.get_neighbor(system, b), Peridynamics.dims(system))
          for b in Peridynamics.each_bond_idx(system, 1)]
    @test sum(Φ[k][1] * ΔX[k][1] for k in eachindex(Φ)) ≈ 1
    # the deformation gradient of a uniform stretch in x is recovered without `NaN`
    storage.position[1, :] .*= 1.01
    Peridynamics.calc_force_density!(dh, 0.0, 0.0)
    F = Peridynamics.get_tensor(storage.defgrad, 1, Peridynamics.dims(system))
    @test !Peridynamics.containsnan(F)
    @test F[1, 1] ≈ 1.01
    @test all(isfinite, storage.b_int)
end

@testitem "export_field: the strain energy density is computed on demand" setup=[Fixtures] begin
    import Peridynamics: export_field

    body = Fixtures.cube(RKCMaterial(); n=4)
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    chunk.storage.position .*= 1.001 # uniform stretch, so the energy is nonzero
    Peridynamics.calc_weights_and_defgrad!(chunk, 0.0, 1e-7)
    Peridynamics.calc_force_density!(chunk, 0.0, 1e-7)
    export_field(Val(:strain_energy_density), chunk.mat, chunk.system, chunk.storage,
                 chunk.paramsetup, 0.0)
    sed = chunk.storage.strain_energy_density
    @test all(>(0), sed)
    @test all(isfinite, sed)
end
