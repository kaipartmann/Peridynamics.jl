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
