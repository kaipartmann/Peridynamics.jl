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
    D_solution = @SMatrix [0.09091326258480677 -0.004549586777144192 0.0;
                           -0.004549586777144192 0.09999590721233696 0.0;
                           0.0 0.0 0.1]
    @test Dᵣ ≈ D_solution atol=10eps()
end
