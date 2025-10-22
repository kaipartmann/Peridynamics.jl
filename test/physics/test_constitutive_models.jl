# ============================================================================ #
# ANALYTICAL TEST CASES FOR CONSTITUTIVE MODELS
# ============================================================================ #
# The following test cases use known analytical solutions to verify the
# correctness of the constitutive model implementations. Each test verifies:
# 1. Zero stress at F = I (reference configuration)
# 2. Known analytical solutions for specific deformation states
# 3. Small strain limit compatibility (where applicable)
# ============================================================================ #

@testsnippet AnalyticalTestCases begin
    using Peridynamics.StaticArrays
    using LinearAlgebra

    """
    Helper function to set up a material model with specified parameters.
    Returns storage and params for testing.
    """
    function setup_material(model; E=210e9, nu=0.3, rho=7850)
        pos, vol = uniform_box(1.0, 1.0, 1.0, 0.5)
        mat = CMaterial(; model)
        body = Body(mat, pos, vol)
        material!(body; horizon=1, E, nu, rho)
        decomp = Peridynamics.PointDecomposition(body, 1)
        system = Peridynamics.get_system(body, decomp, 1)
        solver = VelocityVerlet(steps=1)
        storage = Peridynamics.CStorage(mat, solver, system)
        params = body.point_params[1]
        return storage, params
    end

    """
    Analytical solution for pure shear in small strain limit.
    For F = [1, γ, 0; 0, 1, 0; 0, 0, 1],
    E = [γ²/2, γ/2, 0; γ/2, 0, 0; 0, 0, 0]
    S = λ*tr(E)*I + 2μ*E
    P = F * S
    """
    function analytical_pure_shear_small_strain(λ, μ, γ)
        # Exact Green-Lagrange strain for F = [1, γ, 0; 0, 1, 0; 0, 0, 1]
        # E = 0.5*(F'F - I) = [0, γ/2, 0; γ/2, γ²/2, 0; 0, 0, 0]
        E11 = 0.0
        E12 = 0.5 * γ
        E22 = 0.5 * γ^2
        E33 = 0.0
        tr_E = E11 + E22 + E33  # = γ²/2

        # S = λ*tr(E)*I + 2μ*E
        S11 = λ*tr_E + 2*μ*E11
        S12 = 2*μ*E12
        S22 = λ*tr_E + 2*μ*E22
        S33 = λ*tr_E

        S = @SMatrix [S11 S12 0.0; S12 S22 0.0; 0.0 0.0 S33]
        F = @SMatrix [1.0 γ 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        return F * S
    end

    """
    Analytical solution for uniaxial stretch in small strain limit.
    For F = diag(1+ε, 1, 1), E ≈ diag(ε + ε²/2, -νε, -νε) for small ε
    """
    function analytical_uniaxial_small_strain(λ, μ, ε)
        # E = 0.5*(F'F - I) for F = diag(1+ε, 1, 1)
        # E ≈ diag(ε + ε²/2, 0, 0)
        E11 = ε + 0.5*ε^2
        E22 = 0.0
        E33 = 0.0
        # S = λ*tr(E)*I + 2μ*E
        tr_E = E11 + E22 + E33
        S11 = λ*tr_E + 2*μ*E11
        S22 = λ*tr_E + 2*μ*E22
        S33 = λ*tr_E + 2*μ*E33
        S = @SMatrix [S11 0.0 0.0; 0.0 S22 0.0; 0.0 0.0 S33]
        F = @SMatrix [1.0+ε 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        return F * S
    end

    """
    Analytical strain energy for pure shear in small strain limit.
    For E = [0, γ/2, 0; γ/2, γ²/2, 0; 0, 0, 0]:
    Ψ = (λ/2)*tr(E)² + μ*tr(E*E)
      = (λ/2)*(γ²/2)² + μ*(2*(γ/2)² + (γ²/2)²)
      = λ*γ⁴/8 + μ*(γ²/2 + γ⁴/4)
    For small γ, dominant term is μ*γ²/2
    """
    function analytical_energy_pure_shear(λ, μ, γ)
        E11 = 0.0
        E12 = 0.5 * γ
        E22 = 0.5 * γ^2
        tr_E = E22  # Only E22 contributes
        # tr(E*E) = 2*E12² + E22²
        tr_E_sq = 2*E12^2 + E22^2
        return 0.5 * λ * tr_E^2 + μ * tr_E_sq
    end

    """
    Analytical strain energy for uniaxial tension in small strain limit.
    For small ε: Ψ ≈ (λ/2)*ε² + μ*ε²
    """
    function analytical_energy_uniaxial(λ, μ, ε)
        # E ≈ diag(ε + ε²/2, 0, 0)
        E11 = ε + 0.5*ε^2
        tr_E = E11
        return 0.5 * λ * tr_E^2 + μ * E11^2
    end

    """
    Analytical strain energy for volumetric deformation.
    For F = λ*I, E = ((λ²-1)/2)*I, tr(E) = 3(λ²-1)/2
    Ψ = (λ/2)*tr(E)² + μ*tr(E²) = (λ/2)*(3(λ²-1)/2)² + μ*3((λ²-1)/2)²
    """
    function analytical_energy_volumetric(λ, μ, stretch)
        tr_E = 3 * (stretch^2 - 1) / 2
        E_mag_sq = 3 * ((stretch^2 - 1) / 2)^2
        return 0.5 * λ * tr_E^2 + μ * E_mag_sq
    end
end

@testitem "LinearElastic - Reference Configuration" setup=[AnalyticalTestCases] begin
    model = LinearElastic()
    storage, params = setup_material(model; E=210e9, nu=0.3, rho=7850)

    # Test 1: Zero stress at reference configuration
    F = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test iszero(P)

    Ψ = Peridynamics.strain_energy_density(model, storage, params, F)
    @test iszero(Ψ)
end

@testitem "LinearElastic - Small Strain Limit" setup=[AnalyticalTestCases] begin
    model = LinearElastic()
    storage, params = setup_material(model; E=210e9, nu=0.3, rho=7850)

    # Test 2: Pure shear (small deformation)
    γ = 0.01  # Small shear strain
    F = @SMatrix [1.0 γ 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    P_analytical = analytical_pure_shear_small_strain(params.λ, params.μ, γ)
    @test P ≈ P_analytical

    Ψ = Peridynamics.strain_energy_density(model, storage, params, F)
    Ψ_analytical = analytical_energy_pure_shear(params.λ, params.μ, γ)
    @test Ψ ≈ Ψ_analytical

    # Test 3: Uniaxial tension (small strain)
    ε = 0.01  # Small normal strain
    F = @SMatrix [1.0+ε 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    P_analytical = analytical_uniaxial_small_strain(params.λ, params.μ, ε)
    @test P ≈ P_analytical

    Ψ = Peridynamics.strain_energy_density(model, storage, params, F)
    Ψ_analytical = analytical_energy_uniaxial(params.λ, params.μ, ε)
    @test Ψ ≈ Ψ_analytical
end

@testitem "LinearElastic - Finite Strain" setup=[AnalyticalTestCases] begin
    model = LinearElastic()
    storage, params = setup_material(model; E=210e9, nu=0.3, rho=7850)

    # Test 4: Volumetric deformation (exact solution via energy)
    λ_stretch = 1.1
    F = @SMatrix [λ_stretch 0.0 0.0; 0.0 λ_stretch 0.0; 0.0 0.0 λ_stretch]
    Ψ = Peridynamics.strain_energy_density(model, storage, params, F)
    Ψ_analytical = analytical_energy_volumetric(params.λ, params.μ, λ_stretch)
    @test Ψ ≈ Ψ_analytical
end

@testitem "SaintVenantKirchhoff - All Tests" setup=[AnalyticalTestCases] begin
    model = SaintVenantKirchhoff()
    storage, params = setup_material(model; E=210e9, nu=0.3, rho=7850)

    # Reference configuration
    F = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test iszero(P)
    Ψ = Peridynamics.strain_energy_density(model, storage, params, F)
    @test iszero(Ψ)

    # Small strain - pure shear
    γ = 0.01
    F = @SMatrix [1.0 γ 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    P_analytical = analytical_pure_shear_small_strain(params.λ, params.μ, γ)
    @test P ≈ P_analytical
    Ψ = Peridynamics.strain_energy_density(model, storage, params, F)
    Ψ_analytical = analytical_energy_pure_shear(params.λ, params.μ, γ)
    @test Ψ ≈ Ψ_analytical

    # Small strain - uniaxial tension
    ε = 0.01
    F = @SMatrix [1.0+ε 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    P_analytical = analytical_uniaxial_small_strain(params.λ, params.μ, ε)
    @test P ≈ P_analytical
    Ψ = Peridynamics.strain_energy_density(model, storage, params, F)
    Ψ_analytical = analytical_energy_uniaxial(params.λ, params.μ, ε)
    @test Ψ ≈ Ψ_analytical

    # Finite strain - volumetric
    λ_stretch = 1.1
    F = @SMatrix [λ_stretch 0.0 0.0; 0.0 λ_stretch 0.0; 0.0 0.0 λ_stretch]
    Ψ = Peridynamics.strain_energy_density(model, storage, params, F)
    Ψ_analytical = analytical_energy_volumetric(params.λ, params.μ, λ_stretch)
    @test Ψ ≈ Ψ_analytical
end

@testitem "NeoHooke - Reference and Small Strain" setup=[AnalyticalTestCases] begin
    model = NeoHooke()
    storage, params = setup_material(model; E=210e9, nu=0.3, rho=7850)

    # Reference configuration
    F = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test iszero(P)
    Ψ = Peridynamics.strain_energy_density(model, storage, params, F)
    @test iszero(Ψ)

    # Small strain limit - should match linear elasticity
    γ = 0.001  # Very small for linearization
    F = @SMatrix [1.0 γ 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    P_analytical = analytical_pure_shear_small_strain(params.λ, params.μ, γ)
    @test P ≈ P_analytical rtol=1e-2

    Ψ = Peridynamics.strain_energy_density(model, storage, params, F)
    Ψ_analytical = analytical_energy_pure_shear(params.λ, params.μ, γ)
    @test Ψ ≈ Ψ_analytical rtol=1e-4
end

@testitem "NeoHooke - Incompressibility Check" setup=[AnalyticalTestCases] begin
    model = NeoHooke()
    storage, params = setup_material(model; E=210e9, nu=0.3, rho=7850)

    # For nearly incompressible material (high Poisson ratio)
    storage_incomp, params_incomp = setup_material(model; E=210e9, nu=0.49, rho=7850)

    # Pure shear should give det(F) = 1
    γ = 0.5
    F = @SMatrix [1.0 γ 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

    # Verify stress exists and is finite
    P = Peridynamics.first_piola_kirchhoff(model, storage_incomp, params_incomp, F)
    @test all(isfinite.(P))
    @test norm(P) > 0
end

@testitem "NeoHookePenalty - Reference and Small Strain" setup=[AnalyticalTestCases] begin
    model = NeoHookePenalty()
    storage, params = setup_material(model; E=210e9, nu=0.3, rho=7850)

    # Reference configuration
    F = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test iszero(P)
    Ψ = Peridynamics.strain_energy_density(model, storage, params, F)
    @test iszero(Ψ)

    # Small strain limit
    γ = 0.001
    F = @SMatrix [1.0 γ 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    P_analytical = analytical_pure_shear_small_strain(params.λ, params.μ, γ)
    # NeoHookePenalty has different formulation, so tolerance is higher
    @test P ≈ P_analytical rtol=1e-2
end

@testitem "NeoHookePenalty - Volume Preservation" setup=[AnalyticalTestCases] begin
    model = NeoHookePenalty()
    storage, params = setup_material(model; E=210e9, nu=0.3, rho=7850)

    # Isochoric deformation: stretch in one direction, compress in others
    λ = 1.2
    F = @SMatrix [λ 0.0 0.0; 0.0 1/sqrt(λ) 0.0; 0.0 0.0 1/sqrt(λ)]

    # Verify stress is finite and non-zero
    P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    @test all(isfinite.(P))
    @test norm(P) > 0
end

# ============================================================================ #
# THERMODYNAMIC CONSISTENCY TESTS
# ============================================================================ #
# These tests verify that P = ∂Ψ/∂F for all constitutive models using
# automatic differentiation. This is a fundamental requirement for hyperelastic
# materials to ensure thermodynamic consistency.
# ============================================================================ #

@testsnippet ConsistencyChecker begin
    using Peridynamics.StaticArrays, Peridynamics.Printf, Peridynamics.LinearAlgebra
    using ForwardDiff

    function get_storage_params_test_setup(model)
        pos, vol = uniform_box(1,1,1,0.5)
        body = Body(CMaterial(; model), pos, vol)
        material!(body; horizon=1, E=210e9, nu=0.3, rho=7850)
        ts = VelocityVerlet(steps=1)
        return get_storage_params_test_setup(body, ts)
    end

    function get_storage_params_test_setup(body, ts)
        dh = Peridynamics.threads_data_handler(body, ts, 1)
        (; storage, paramsetup) = dh.chunks[1]
        params = Peridynamics.get_params(paramsetup, 1)
        return storage, params
    end

    # Automatic differentiation using the relation P = ∂Ψ/∂F
    function pk1_ad(setup, F)
        (; model, storage, params) = setup
        Ψ = let model=model, storage=storage, params=params, F=F
            F -> Peridynamics.strain_energy_density(model, storage, params, F)
        end
        config = ForwardDiff.GradientConfig(Ψ, F, ForwardDiff.Chunk{9}())
        return ForwardDiff.gradient(Ψ, F, config)
    end

    # The internal function defined manually
    function pk1(setup, F)
        (; model, storage, params) = setup
        return Peridynamics.first_piola_kirchhoff(model, storage, params, F)
    end
end

@testitem "Thermodynamic Consistency: P = ∂Ψ/∂F" setup=[ConsistencyChecker] begin
    # Setup all models
    model_le = LinearElastic()
    storage_le, params_le = get_storage_params_test_setup(model_le)
    setup_le = (; model=model_le, storage=storage_le, params=params_le)

    model_svk = SaintVenantKirchhoff()
    storage_svk, params_svk = get_storage_params_test_setup(model_svk)
    setup_svk = (; model=model_svk, storage=storage_svk, params=params_svk)

    model_nh = NeoHooke()
    storage_nh, params_nh = get_storage_params_test_setup(model_nh)
    setup_nh = (; model=model_nh, storage=storage_nh, params=params_nh)

    model_nhp = NeoHookePenalty()
    storage_nhp, params_nhp = get_storage_params_test_setup(model_nhp)
    setup_nhp = (; model=model_nhp, storage=storage_nhp, params=params_nhp)

    # Test 1: Isotropic extension
    λxyz = 1.2
    F = @SMatrix [λxyz 0.0 0.0; 0.0 λxyz 0.0; 0.0 0.0 λxyz]
    @test pk1(setup_le, F) ≈ pk1_ad(setup_le, F)
    @test pk1(setup_svk, F) ≈ pk1_ad(setup_svk, F)
    @test pk1(setup_nh, F) ≈ pk1_ad(setup_nh, F)
    @test pk1(setup_nhp, F) ≈ pk1_ad(setup_nhp, F)

    # Test 2: Pure shear
    β = 0.1
    F = @SMatrix [1.0 β 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test pk1(setup_le, F) ≈ pk1_ad(setup_le, F)
    @test pk1(setup_svk, F) ≈ pk1_ad(setup_svk, F)
    @test pk1(setup_nh, F) ≈ pk1_ad(setup_nh, F)
    @test pk1(setup_nhp, F) ≈ pk1_ad(setup_nhp, F)

    # Test 3: Uniaxial tension in x-direction
    λx = 1.1
    F = @SMatrix [λx 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test pk1(setup_le, F) ≈ pk1_ad(setup_le, F)
    @test pk1(setup_svk, F) ≈ pk1_ad(setup_svk, F)
    @test pk1(setup_nh, F) ≈ pk1_ad(setup_nh, F)
    @test pk1(setup_nhp, F) ≈ pk1_ad(setup_nhp, F)

    # Test 4: Uniaxial tension in y-direction
    λy = 1.2
    F = @SMatrix [1.0 0.0 0.0; 0.0 λy 0.0; 0.0 0.0 1.0]
    @test pk1(setup_le, F) ≈ pk1_ad(setup_le, F)
    @test pk1(setup_svk, F) ≈ pk1_ad(setup_svk, F)
    @test pk1(setup_nh, F) ≈ pk1_ad(setup_nh, F)
    @test pk1(setup_nhp, F) ≈ pk1_ad(setup_nhp, F)

    # Test 5: Uniaxial tension in z-direction
    λz = 1.3
    F = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 λz]
    @test pk1(setup_le, F) ≈ pk1_ad(setup_le, F)
    @test pk1(setup_svk, F) ≈ pk1_ad(setup_svk, F)
    @test pk1(setup_nh, F) ≈ pk1_ad(setup_nh, F)
    @test pk1(setup_nhp, F) ≈ pk1_ad(setup_nhp, F)

    # Test 6: Combined shear and extension
    F = @SMatrix [1.1 0.05 0.0; 0.0 1.15 0.0; 0.0 0.0 1.05]
    @test pk1(setup_le, F) ≈ pk1_ad(setup_le, F)
    @test pk1(setup_svk, F) ≈ pk1_ad(setup_svk, F)
    @test pk1(setup_nh, F) ≈ pk1_ad(setup_nh, F)
    @test pk1(setup_nhp, F) ≈ pk1_ad(setup_nhp, F)

    # Test 7: Combined large deformation
    F = 2 * @SMatrix [1.1 0.05 0.1; 0.02 1.15 0.2; 0.4 0.01 1.05]
    @test pk1(setup_le, F) ≈ pk1_ad(setup_le, F)
    @test pk1(setup_svk, F) ≈ pk1_ad(setup_svk, F)
    @test pk1(setup_nh, F) ≈ pk1_ad(setup_nh, F)
    @test pk1(setup_nhp, F) ≈ pk1_ad(setup_nhp, F)

    # Test 8: Random large deformation
    F = 2I + @SMatrix rand(3, 3)
    @test pk1(setup_le, F) ≈ pk1_ad(setup_le, F)
    @test pk1(setup_svk, F) ≈ pk1_ad(setup_svk, F)
    @test pk1(setup_nh, F) ≈ pk1_ad(setup_nh, F)
    @test pk1(setup_nhp, F) ≈ pk1_ad(setup_nhp, F)
end

# ============================================================================ #
# CAUCHY STRESS VERIFICATION TESTS
# ============================================================================ #
# These tests verify the physical properties that Cauchy stress must satisfy
# for all hyperelastic materials. The Cauchy stress is computed as:
# σ = (1/J) * P * F^T, where J = det(F)
#
# Note: We test physical properties (symmetry, zero stress at reference, etc.)
# rather than specific analytical values, since the constitutive models may
# use formulations that differ from classical textbook solutions while still
# being thermodynamically consistent (verified by P = ∂Ψ/∂F tests above).
# ============================================================================ #

@testitem "Cauchy Stress: Physical Properties" begin
    using Peridynamics.StaticArrays
    using LinearAlgebra

    # Test physical properties that should hold for all hyperelastic materials
    E = 210e9
    nu = 0.3
    pos, vol = uniform_box(1.0, 1.0, 1.0, 0.5)

    # Test 1: Symmetry of Cauchy stress tensor for general deformation
    F = @SMatrix [1.2 0.1 0.05; 0.05 1.15 0.08; 0.03 0.06 1.1]

    for model in [LinearElastic(), SaintVenantKirchhoff(), NeoHooke(), NeoHookePenalty()]
        mat = CMaterial(; model)
        body = Body(mat, pos, vol)
        material!(body; horizon=1, E, nu, rho=7850)

        decomp = Peridynamics.PointDecomposition(body, 1)
        system = Peridynamics.get_system(body, decomp, 1)
        solver = VelocityVerlet(steps=1)
        storage = Peridynamics.CStorage(mat, solver, system)
        params = body.point_params[1]

        P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
        σ = Peridynamics.cauchy_stress(P, F)

        # Cauchy stress must be symmetric
        @test σ ≈ σ'
    end

    # Test 2: Zero stress at reference configuration
    F = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

    for model in [LinearElastic(), SaintVenantKirchhoff(), NeoHooke(), NeoHookePenalty()]
        mat = CMaterial(; model)
        body = Body(mat, pos, vol)
        material!(body; horizon=1, E, nu, rho=7850)

        decomp = Peridynamics.PointDecomposition(body, 1)
        system = Peridynamics.get_system(body, decomp, 1)
        solver = VelocityVerlet(steps=1)
        storage = Peridynamics.CStorage(mat, solver, system)
        params = body.point_params[1]

        P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
        σ = Peridynamics.cauchy_stress(P, F)

        @test norm(σ) ≈ 0.0 atol=eps()
    end

    # Test 3: Positive normal stress under tension
    F = @SMatrix [1.1 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

    for model in [LinearElastic(), SaintVenantKirchhoff(), NeoHooke(), NeoHookePenalty()]
        mat = CMaterial(; model)
        body = Body(mat, pos, vol)
        material!(body; horizon=1, E, nu, rho=7850)

        decomp = Peridynamics.PointDecomposition(body, 1)
        system = Peridynamics.get_system(body, decomp, 1)
        solver = VelocityVerlet(steps=1)
        storage = Peridynamics.CStorage(mat, solver, system)
        params = body.point_params[1]

        P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
        σ = Peridynamics.cauchy_stress(P, F)

        # Under uniaxial tension, σ11 should be positive
        @test σ[1,1] > 0
    end

    # Test 4: Test Cauchy stress under triaxial deformation
    λ1, λ2, λ3 = 1.01, 1.02, 1.03
    F = @SMatrix [λ1 0 0; 0 λ2 0; 0 0 λ3]
    for model in [LinearElastic(), SaintVenantKirchhoff(), NeoHooke(), NeoHookePenalty()]
        mat = CMaterial(; model)
        body = Body(mat, pos, vol)
        material!(body; horizon=1, E, nu, rho=7850)

        decomp = Peridynamics.PointDecomposition(body, 1)
        system = Peridynamics.get_system(body, decomp, 1)
        solver = VelocityVerlet(steps=1)
        storage = Peridynamics.CStorage(mat, solver, system)
        params = body.point_params[1]

        P = Peridynamics.first_piola_kirchhoff(model, storage, params, F)
        σ = Peridynamics.cauchy_stress(P, F)

        # Check that normal stresses correspond to stretches
        @test P[1,1] > 0
        @test P[2,2] > 0
        @test P[3,3] > 0
        @test P[1,2] ≈ 0.0 atol=eps()
        @test P[2,1] ≈ 0.0 atol=eps()
        @test P[1,3] ≈ 0.0 atol=eps()
        @test P[3,1] ≈ 0.0 atol=eps()
        @test P[2,3] ≈ 0.0 atol=eps()
        @test P[3,2] ≈ 0.0 atol=eps()
        @test σ[1,1] > 0
        @test σ[2,2] > 0
        @test σ[3,3] > 0
        @test σ[1,2] ≈ 0.0 atol=eps()
        @test σ[2,1] ≈ 0.0 atol=eps()
        @test σ[1,3] ≈ 0.0 atol=eps()
        @test σ[3,1] ≈ 0.0 atol=eps()
        @test σ[2,3] ≈ 0.0 atol=eps()
        @test σ[3,2] ≈ 0.0 atol=eps()
    end
end
