# The exported Cauchy, von Mises and hydrostatic stress of the correspondence materials under
# homogeneous deformations, checked against the small strain closed forms. One item looping
# `Fixtures.MODEL_MATERIALS`; adding a (material, model) pair there covers it here.

@testmodule StressExport begin
    using Peridynamics, Test
    using Peridynamics.StaticArrays, Peridynamics.LinearAlgebra

    # the exported field `field` of `body` after the homogeneous deformation `F`
    function export_field(body, F, field)
        ts = VelocityVerlet(steps=1)
        dh = Peridynamics.threads_data_handler(body, ts, 1)
        Peridynamics.init_time_solver!(ts, dh)
        Peridynamics.initialize!(dh, ts)
        (; mat, system, storage, paramsetup) = dh.chunks[1]
        for i in Peridynamics.each_point_idx(system)
            Xi = Peridynamics.get_vector(system.position, i)
            xi = F * Xi
            Peridynamics.update_vector!(storage.position, i, xi)
            vi = (xi - Xi) / ts.Δt
            Peridynamics.update_vector!(storage.velocity, i, vi)
            Peridynamics.update_vector!(storage.velocity_half, i, vi)
        end
        Peridynamics.calc_force_density!(dh, ts.Δt, ts.Δt)
        return Peridynamics.export_field(Val(field), mat, system, storage, paramsetup, 0.0)
    end

    # mean of the `3×3` tensors stored column-wise in the `9×n` matrix `qty`
    function mean_tensor(qty)
        s = zero(SMatrix{3,3,eltype(qty),9})
        for i in axes(qty, 2)
            s += Peridynamics.get_tensor(qty, i)
        end
        return s / size(qty, 2)
    end

    const OFF_DIAGONAL = ((1, 2), (2, 3), (3, 1), (2, 1), (3, 2), (1, 3))

    # check the stress exports of `body` under a homogeneous isotropic extension, a pure shear
    # and three uniaxial extensions against the small strain continuum with the point parameters
    # of `body`; `rtol` is the relative tolerance on the normal and hydrostatic stresses (the
    # Neo-Hooke models deviate from the linear closed forms by their nonlinearity)
    function check(body; rtol=1e-3)
        p = body.point_params[1]
        λ = 1.001 # uniform extension
        β = 0.001 # shear parameter

        # homogeneous isotropic extension: σ = 3K(λ-1) I, no deviatoric part
        F = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
        σ = mean_tensor(export_field(body, F, :cauchy_stress))
        σ_ii = 3 * p.K * (λ - 1)
        @test all(isapprox(σ[i, i], σ_ii; rtol) for i in 1:3)
        @test all(isapprox(σ[i, j], 0.0; atol=1e-5) for (i, j) in OFF_DIAGONAL)
        @test all(export_field(body, F, :von_mises_stress) .< 1e-7 * σ_ii)
        @test all(isapprox.(export_field(body, F, :hydrostatic_stress), σ_ii; rtol))

        # pure shear: σ₁₂ = μβ, no normal and no hydrostatic stress
        F = @SMatrix [1 β 0; 0 1 0; 0 0 1]
        σ = mean_tensor(export_field(body, F, :cauchy_stress))
        σ_12 = p.μ * β
        @test σ[1, 2] ≈ σ_12 rtol=1e-3
        @test σ[2, 1] ≈ σ_12 rtol=1e-3
        @test all(σ[i, i] < 0.005 * σ_12 for i in 1:3)
        @test all(isapprox(σ[i, j], 0.0; atol=1e-5) for (i, j) in ((2, 3), (3, 1), (3, 2), (1, 3)))
        @test all(isapprox.(export_field(body, F, :von_mises_stress), p.μ * β * √3; rtol=1e-3))
        @test all(isapprox.(export_field(body, F, :hydrostatic_stress), 0; atol=0.005σ_12))

        # uniaxial extension in each direction: σ_dd = (λ+2μ)ε, the others λε
        for d in 1:3
            F = SMatrix{3,3}(i == j ? (i == d ? λ : 1.0) : 0.0 for i in 1:3, j in 1:3)
            σ = mean_tensor(export_field(body, F, :cauchy_stress))
            σ_dd = (p.λ + 2 * p.μ) * (λ - 1)
            σ_other = p.λ * (λ - 1)
            @test all(isapprox(σ[i, i], i == d ? σ_dd : σ_other; rtol=1e-2) for i in 1:3)
            @test all(isapprox(σ[i, j], 0.0; atol=1e-5) for (i, j) in OFF_DIAGONAL)
            @test all(isapprox.(export_field(body, F, :von_mises_stress), abs(σ_dd - σ_other);
                                rtol=1e-2))
            @test all(isapprox.(export_field(body, F, :hydrostatic_stress),
                                (σ_dd + 2σ_other) / 3; rtol))
        end
        return nothing
    end
end

@testitem "stress export: correspondence materials" tags=[:slow] setup=[Fixtures, StressExport] begin
    # `:slow`: 12 material/model combinations, each compiled once; runs in the `extras` CI job
    Δx = 0.25
    for (name, mat) in Fixtures.MODEL_MATERIALS
        # the bond-associated formulation does not export the Cauchy stress
        mat isa BACMaterial && continue
        @testset "$name" begin
            pos, vol = uniform_box(1, 1, 1, Δx)
            body = Body(mat, pos, vol)
            material!(body; horizon=3.01Δx, rho=8000, E=210e9, nu=0.25)
            rtol = mat.constitutive_model isa Union{NeoHooke,NeoHookePenalty} ? 1e-2 : 1e-3
            StressExport.check(body; rtol)
        end
    end
end
