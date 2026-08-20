# The zero-energy mode stabilizations of the correspondence formulations
# (`src/discretization/zem_stabilization.jl`).

@testitem "ZEMSilling and ZEMWan: construction" begin
    @test ZEMSilling().Cs == 0.5
    @test ZEMSilling(Cs=2).Cs == 2.0
    @test ZEMWan() isa Peridynamics.AbstractZEMStabilization
    # the stabilization is the `correction` of a correspondence material
    mat = CMaterial(zem=ZEMWan())
    @test Peridynamics.get_correction(mat, 1, 1, 1) === mat.zem
    @test CMaterial().zem isa ZEMSilling
end

@testitem "calc_zem_stiffness_tensor: contraction of the stiffness with the shape tensor" begin
    # C₁ᵢⱼ = Σₖₗ Cᵢⱼₖₗ K⁻¹ₖₗ; for the isotropic stiffness C = λ δᵢⱼδₖₗ + μ(δᵢₖδⱼₗ + δᵢₗδⱼₖ)
    # this is λ tr(K⁻¹) I + 2μ K⁻¹, and a rotation of the stiffness by the identity changes
    # nothing
    using Peridynamics.StaticArrays, Peridynamics.LinearAlgebra
    λ, μ = 2.0, 3.0
    C = zeros(3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        C[i, j, k, l] = λ * (i == j) * (k == l) + μ * ((i == k) * (j == l) + (i == l) * (j == k))
    end
    Kinv = @SMatrix [2.0 0.1 0.0; 0.1 1.0 0.3; 0.0 0.3 0.5]
    C1 = Peridynamics.calc_zem_stiffness_tensor(C, Kinv)
    @test C1 ≈ λ * tr(Kinv) * I + 2μ * Kinv
    C_rotated = zeros(3, 3, 3, 3)
    R = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test Peridynamics.calc_rotated_zem_stiffness_tensor!(C_rotated, C, Kinv, R) ≈ C1
    @test C_rotated ≈ C
    # an isotropic stiffness is invariant under any rotation
    θ = 0.3
    R = @SMatrix [cos(θ) -sin(θ) 0.0; sin(θ) cos(θ) 0.0; 0.0 0.0 1.0]
    @test Peridynamics.calc_rotated_zem_stiffness_tensor!(C_rotated, C, Kinv, R) ≈ C1
end
