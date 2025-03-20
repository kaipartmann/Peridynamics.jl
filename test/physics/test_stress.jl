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
