# `RKCRMaterial`: the reproducing kernel correspondence formulation with stress rotation of
# `src/physics/rk_correspondence_rotated.jl`.

@testitem "RKCRMaterial: constructor and supported models" begin
    # Test default constructor
    mat1 = RKCRMaterial()
    @test mat1.kernel == const_one_kernel
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
    @test_throws ArgumentError RKCRMaterial(lambda = -0.5)
    @test_throws ArgumentError RKCRMaterial(beta = -0.5)
end

@testitem "RKCRMaterial: deformation gradient and its rate" begin
    using Peridynamics.StaticArrays, Peridynamics.LinearAlgebra
    pos, vol = uniform_box(1, 1, 1, 0.25)
    body = Body(RKCRMaterial(), pos, vol)
    material!(body; horizon=0.76, rho=1, E=210e9, nu=0.25, Gc=1.0)
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    (; storage) = dh.chunks[1]
    (; position, defgrad, defgrad_dot, velocity_half) = storage
    # no displacement: the identity
    Peridynamics.calc_force_density!(dh, 0.0, 0.0)
    @test all(isapprox(Peridynamics.get_tensor(defgrad, i), I; atol=1e-12) for i in eachindex(vol))
    # a small uniform stretch in x
    F_a = @SMatrix [1.00001 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    for i in eachindex(vol)
        Peridynamics.update_vector!(position, i, F_a * Peridynamics.get_vector(position, i))
    end
    Peridynamics.calc_force_density!(dh, 0.0, 0.0)
    @test all(isapprox(Peridynamics.get_tensor(defgrad, i), F_a; atol=1e-5) for i in eachindex(vol))
    # a velocity field v = 0.1 x e₁ has the rate Ḟ = 0.1 F₁₁ e₁ ⊗ e₁, nothing else
    velocity_half .= 0.0
    for i in eachindex(vol)
        velocity_half[1, i] = position[1, i] * 0.1
    end
    Peridynamics.calc_force_density!(dh, 0.0, 0.0)
    @test all(eachindex(vol)) do i
        F_dot = Peridynamics.get_tensor(defgrad_dot, i)
        F_dot[1, 1] > 0 && abs(F_dot[2, 2]) < 1e-14 && abs(F_dot[3, 3]) < 1e-14
    end
end
