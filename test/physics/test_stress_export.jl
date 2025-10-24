@testsnippet StressExport begin
    using Peridynamics.StaticArrays, Peridynamics.LinearAlgebra
    function get_export_field(body, ts, F, field)
        dh = Peridynamics.threads_data_handler(body, ts, 1)
        Peridynamics.init_time_solver!(ts, dh)
        Peridynamics.initialize!(dh, ts)
        (; mat, system, storage, paramsetup) = dh.chunks[1]
        # apply deformation gradient F
        for i in Peridynamics.each_point_idx(system)
            Xi = Peridynamics.get_vector(system.position, i)
            xi = F * Xi
            Peridynamics.update_vector!(storage.position, i, xi)
            vi = (xi - Xi) / ts.Δt
            Peridynamics.update_vector!(storage.velocity, i, vi)
            Peridynamics.update_vector!(storage.velocity_half, i, vi)
        end
        Peridynamics.calc_force_density!(dh, ts.Δt, ts.Δt)
        qty = Peridynamics.export_field(Val(field), mat, system, storage, paramsetup, 0.0)
        return qty
    end
    function mean_tensor(qty)
        s = zero(SMatrix{3,3,eltype(qty),9})
        N = size(qty, 2)
        for i in axes(qty, 2)
            s += Peridynamics.get_tensor(qty, i)
        end
        return s / N
    end
    mean(x) = sum(x) / length(x)
end

@testitem "Exported stress CMaterial SaintVenantKirchhoff" setup=[StressExport] begin
    Δx = 0.25
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = CMaterial(; model=SaintVenantKirchhoff())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.001 # shear parameter

    #===== homogeneous isotropic extension =====#
    F = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_ii_exact = 3 * params.K * (λ - 1)
    # test normal stresses
    for i in 1:3
        @test σ_pd_mean[i,i] ≈ σ_ii_exact rtol=1e-3
    end
    # test shear stresses
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    # should be nearly zero for isotropic extension, so we test that it is
    # at least 7 orders of magnitude smaller than normal stress
    @test all(σ_pd_vm_vec .< 1e-7 * σ_ii_exact)
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    @test all(isapprox.(σ_pd_hydro_vec, σ_ii_exact; rtol=1e-3))

    #===== pure shear =====#
    F = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_12_21_exact = params.μ * β
    # test affected shear components
    for (i, j) in ((1,2), (2,1))
        @test σ_pd_mean[i,j] ≈ σ_12_21_exact rtol=1e-3
    end
    # test normal stresses, should be at least small compared to shear stress
    for i in 1:3
        @test σ_pd_mean[i,i] < 0.005 * σ_12_21_exact
    end
    # test other shear stresses, should be zero
    for (i, j) in ((2,3), (3,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = params.μ * β * √3
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-3))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    # should be nearly zero for pure shear, therefore at least small
    # compared to shear stress
    @test all(isapprox.(σ_pd_hydro_vec, 0; atol=0.005σ_12_21_exact))

    #===== uniform deformation in x =====#
    F = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_11_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_22_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_33_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_22_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_11_exact - σ_22_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_11_exact + 2 * σ_22_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in y =====#
    F = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_22_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_33_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_11_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_22_exact - σ_11_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_22_exact + 2 * σ_11_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in z =====#
    F = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_33_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_22_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_33_exact - σ_11_22_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_33_exact + 2 * σ_11_22_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))
end

@testitem "Exported stress CMaterial LinearElastic" setup=[StressExport] begin
    Δx = 0.25
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = CMaterial(; model=LinearElastic())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.001 # shear parameter

    #===== homogeneous isotropic extension =====#
    F = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_ii_exact = 3 * params.K * (λ - 1)
    # test normal stresses
    for i in 1:3
        @test σ_pd_mean[i,i] ≈ σ_ii_exact rtol=1e-3
    end
    # test shear stresses
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    # should be nearly zero for isotropic extension, so we test that it is
    # at least 7 orders of magnitude smaller than normal stress
    @test all(σ_pd_vm_vec .< 1e-7 * σ_ii_exact)
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    @test all(isapprox.(σ_pd_hydro_vec, σ_ii_exact; rtol=1e-3))

    #===== pure shear =====#
    F = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_12_21_exact = params.μ * β
    # test affected shear components
    for (i, j) in ((1,2), (2,1))
        @test σ_pd_mean[i,j] ≈ σ_12_21_exact rtol=1e-3
    end
    # test normal stresses, should be at least small compared to shear stress
    for i in 1:3
        @test σ_pd_mean[i,i] < 0.005 * σ_12_21_exact
    end
    # test other shear stresses, should be zero
    for (i, j) in ((2,3), (3,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = params.μ * β * √3
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-3))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    # should be nearly zero for pure shear, therefore at least small
    # compared to shear stress
    @test all(isapprox.(σ_pd_hydro_vec, 0; atol=0.005σ_12_21_exact))

    #===== uniform deformation in x =====#
    F = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_11_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_22_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_33_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_22_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_11_exact - σ_22_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_11_exact + 2 * σ_22_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in y =====#
    F = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_22_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_33_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_11_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_22_exact - σ_11_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_22_exact + 2 * σ_11_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in z =====#
    F = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_33_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_22_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_33_exact - σ_11_22_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_33_exact + 2 * σ_11_22_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))
end

@testitem "Exported stress CMaterial NeoHooke" setup=[StressExport] begin
    Δx = 0.25
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = CMaterial(; model=NeoHooke())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.001 # shear parameter

    #===== homogeneous isotropic extension =====#
    F = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_ii_exact = 3 * params.K * (λ - 1)
    # test normal stresses
    for i in 1:3
        @test σ_pd_mean[i,i] ≈ σ_ii_exact rtol=1e-2
    end
    # test shear stresses
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    # should be nearly zero for isotropic extension, so we test that it is
    # at least 7 orders of magnitude smaller than normal stress
    @test all(σ_pd_vm_vec .< 1e-7 * σ_ii_exact)
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    @test all(isapprox.(σ_pd_hydro_vec, σ_ii_exact; rtol=1e-2))

    #===== pure shear =====#
    F = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_12_21_exact = params.μ * β
    # test affected shear components
    for (i, j) in ((1,2), (2,1))
        @test σ_pd_mean[i,j] ≈ σ_12_21_exact rtol=1e-3
    end
    # test normal stresses, should be at least small compared to shear stress
    for i in 1:3
        @test σ_pd_mean[i,i] < 0.005 * σ_12_21_exact
    end
    # test other shear stresses, should be zero
    for (i, j) in ((2,3), (3,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = params.μ * β * √3
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-3))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    # should be nearly zero for pure shear, therefore at least small
    # compared to shear stress
    @test all(isapprox.(σ_pd_hydro_vec, 0; atol=0.005σ_12_21_exact))

    #===== uniform deformation in x =====#
    F = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_11_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_22_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_33_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_22_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_11_exact - σ_22_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_11_exact + 2 * σ_22_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-2))

    #===== uniform deformation in y =====#
    F = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_22_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_33_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_11_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_22_exact - σ_11_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_22_exact + 2 * σ_11_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-2))

    #===== uniform deformation in z =====#
    F = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_33_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_22_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_33_exact - σ_11_22_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_33_exact + 2 * σ_11_22_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-2))
end

@testitem "Exported stress CMaterial NeoHookePenalty" setup=[StressExport] begin
    Δx = 0.25
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = CMaterial(; model=NeoHookePenalty())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.001 # shear parameter

    #===== homogeneous isotropic extension =====#
    F = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_ii_exact = 3 * params.K * (λ - 1)
    # test normal stresses
    for i in 1:3
        @test σ_pd_mean[i,i] ≈ σ_ii_exact rtol=1e-2
    end
    # test shear stresses
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    # should be nearly zero for isotropic extension, so we test that it is
    # at least 7 orders of magnitude smaller than normal stress
    @test all(σ_pd_vm_vec .< 1e-7 * σ_ii_exact)
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    @test all(isapprox.(σ_pd_hydro_vec, σ_ii_exact; rtol=1e-2))

    #===== pure shear =====#
    F = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_12_21_exact = params.μ * β
    # test affected shear components
    for (i, j) in ((1,2), (2,1))
        @test σ_pd_mean[i,j] ≈ σ_12_21_exact rtol=1e-3
    end
    # test normal stresses, should be at least small compared to shear stress
    for i in 1:3
        @test σ_pd_mean[i,i] < 0.005 * σ_12_21_exact
    end
    # test other shear stresses, should be zero
    for (i, j) in ((2,3), (3,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = params.μ * β * √3
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-3))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    # should be nearly zero for pure shear, therefore at least small
    # compared to shear stress
    @test all(isapprox.(σ_pd_hydro_vec, 0; atol=0.005σ_12_21_exact))

    #===== uniform deformation in x =====#
    F = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_11_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_22_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_33_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_22_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_11_exact - σ_22_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_11_exact + 2 * σ_22_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-2))

    #===== uniform deformation in y =====#
    F = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_22_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_33_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_11_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_22_exact - σ_11_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_22_exact + 2 * σ_11_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-2))

    #===== uniform deformation in z =====#
    F = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_33_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_22_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_33_exact - σ_11_22_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_33_exact + 2 * σ_11_22_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-2))
end


@testitem "Exported stress CRMaterial SaintVenantKirchhoff" setup=[StressExport] begin
    Δx = 0.25
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = CRMaterial(; model=SaintVenantKirchhoff())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.001 # shear parameter

    #===== homogeneous isotropic extension =====#
    F = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_ii_exact = 3 * params.K * (λ - 1)
    # test normal stresses
    for i in 1:3
        @test σ_pd_mean[i,i] ≈ σ_ii_exact rtol=1e-3
    end
    # test shear stresses
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    # should be nearly zero for isotropic extension, so we test that it is
    # at least 7 orders of magnitude smaller than normal stress
    @test all(σ_pd_vm_vec .< 1e-7 * σ_ii_exact)
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    @test all(isapprox.(σ_pd_hydro_vec, σ_ii_exact; rtol=1e-3))

    #===== pure shear =====#
    F = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_12_21_exact = params.μ * β
    # test affected shear components
    for (i, j) in ((1,2), (2,1))
        @test σ_pd_mean[i,j] ≈ σ_12_21_exact rtol=1e-3
    end
    # test normal stresses, should be at least small compared to shear stress
    for i in 1:3
        @test σ_pd_mean[i,i] < 0.005 * σ_12_21_exact
    end
    # test other shear stresses, should be zero
    for (i, j) in ((2,3), (3,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = params.μ * β * √3
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-3))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    # should be nearly zero for pure shear, therefore at least small
    # compared to shear stress
    @test all(isapprox.(σ_pd_hydro_vec, 0; atol=0.005σ_12_21_exact))

    #===== uniform deformation in x =====#
    F = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_11_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_22_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_33_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_22_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_11_exact - σ_22_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_11_exact + 2 * σ_22_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in y =====#
    F = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_22_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_33_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_11_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_22_exact - σ_11_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_22_exact + 2 * σ_11_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in z =====#
    F = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_33_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_22_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_33_exact - σ_11_22_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_33_exact + 2 * σ_11_22_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))
end

@testitem "Exported stress CRMaterial LinearElastic" setup=[StressExport] begin
    Δx = 0.25
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = CRMaterial(; model=LinearElastic())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.001 # shear parameter

    #===== homogeneous isotropic extension =====#
    F = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_ii_exact = 3 * params.K * (λ - 1)
    # test normal stresses
    for i in 1:3
        @test σ_pd_mean[i,i] ≈ σ_ii_exact rtol=1e-3
    end
    # test shear stresses
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    # should be nearly zero for isotropic extension, so we test that it is
    # at least 7 orders of magnitude smaller than normal stress
    @test all(σ_pd_vm_vec .< 1e-7 * σ_ii_exact)
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    @test all(isapprox.(σ_pd_hydro_vec, σ_ii_exact; rtol=1e-3))

    #===== pure shear =====#
    F = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_12_21_exact = params.μ * β
    # test affected shear components
    for (i, j) in ((1,2), (2,1))
        @test σ_pd_mean[i,j] ≈ σ_12_21_exact rtol=1e-3
    end
    # test normal stresses, should be at least small compared to shear stress
    for i in 1:3
        @test σ_pd_mean[i,i] < 0.005 * σ_12_21_exact
    end
    # test other shear stresses, should be zero
    for (i, j) in ((2,3), (3,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = params.μ * β * √3
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-3))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    # should be nearly zero for pure shear, therefore at least small
    # compared to shear stress
    @test all(isapprox.(σ_pd_hydro_vec, 0; atol=0.005σ_12_21_exact))

    #===== uniform deformation in x =====#
    F = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_11_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_22_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_33_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_22_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_11_exact - σ_22_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_11_exact + 2 * σ_22_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in y =====#
    F = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_22_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_33_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_11_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_22_exact - σ_11_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_22_exact + 2 * σ_11_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in z =====#
    F = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_33_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_22_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_33_exact - σ_11_22_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_33_exact + 2 * σ_11_22_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))
end

#===========================================================================================
    TESTS FOR THE RKCMaterial
===========================================================================================#
@testitem "Exported stress RKCMaterial SaintVenantKirchhoff" setup=[StressExport] begin
    Δx = 0.25
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = RKCMaterial(; model=SaintVenantKirchhoff())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.001 # shear parameter

    #===== homogeneous isotropic extension =====#
    F = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_ii_exact = 3 * params.K * (λ - 1)
    # test normal stresses
    for i in 1:3
        @test σ_pd_mean[i,i] ≈ σ_ii_exact rtol=1e-3
    end
    # test shear stresses
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    # should be nearly zero for isotropic extension, so we test that it is
    # at least 7 orders of magnitude smaller than normal stress
    @test all(σ_pd_vm_vec .< 1e-7 * σ_ii_exact)
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    @test all(isapprox.(σ_pd_hydro_vec, σ_ii_exact; rtol=1e-3))

    #===== pure shear =====#
    F = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_12_21_exact = params.μ * β
    # test affected shear components
    for (i, j) in ((1,2), (2,1))
        @test σ_pd_mean[i,j] ≈ σ_12_21_exact rtol=1e-3
    end
    # test normal stresses, should be at least small compared to shear stress
    for i in 1:3
        @test σ_pd_mean[i,i] < 0.005 * σ_12_21_exact
    end
    # test other shear stresses, should be zero
    for (i, j) in ((2,3), (3,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = params.μ * β * √3
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-3))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    # should be nearly zero for pure shear, therefore at least small
    # compared to shear stress
    @test all(isapprox.(σ_pd_hydro_vec, 0; atol=0.005σ_12_21_exact))

    #===== uniform deformation in x =====#
    F = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_11_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_22_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_33_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_22_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_11_exact - σ_22_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_11_exact + 2 * σ_22_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in y =====#
    F = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_22_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_33_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_11_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_22_exact - σ_11_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_22_exact + 2 * σ_11_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in z =====#
    F = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_33_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_22_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_33_exact - σ_11_22_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_33_exact + 2 * σ_11_22_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))
end

@testitem "Exported stress RKCMaterial LinearElastic" setup=[StressExport] begin
    Δx = 0.25
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = RKCMaterial(; model=LinearElastic())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.001 # shear parameter

    #===== homogeneous isotropic extension =====#
    F = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_ii_exact = 3 * params.K * (λ - 1)
    # test normal stresses
    for i in 1:3
        @test σ_pd_mean[i,i] ≈ σ_ii_exact rtol=1e-3
    end
    # test shear stresses
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    # should be nearly zero for isotropic extension, so we test that it is
    # at least 7 orders of magnitude smaller than normal stress
    @test all(σ_pd_vm_vec .< 1e-7 * σ_ii_exact)
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    @test all(isapprox.(σ_pd_hydro_vec, σ_ii_exact; rtol=1e-3))

    #===== pure shear =====#
    F = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_12_21_exact = params.μ * β
    # test affected shear components
    for (i, j) in ((1,2), (2,1))
        @test σ_pd_mean[i,j] ≈ σ_12_21_exact rtol=1e-3
    end
    # test normal stresses, should be at least small compared to shear stress
    for i in 1:3
        @test σ_pd_mean[i,i] < 0.005 * σ_12_21_exact
    end
    # test other shear stresses, should be zero
    for (i, j) in ((2,3), (3,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = params.μ * β * √3
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-3))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    # should be nearly zero for pure shear, therefore at least small
    # compared to shear stress
    @test all(isapprox.(σ_pd_hydro_vec, 0; atol=0.005σ_12_21_exact))

    #===== uniform deformation in x =====#
    F = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_11_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_22_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_33_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_22_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_11_exact - σ_22_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_11_exact + 2 * σ_22_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in y =====#
    F = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_22_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_33_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_11_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_22_exact - σ_11_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_22_exact + 2 * σ_11_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in z =====#
    F = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_33_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_22_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_33_exact - σ_11_22_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_33_exact + 2 * σ_11_22_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))
end

@testitem "Exported stress RKCMaterial NeoHooke" setup=[StressExport] begin
    Δx = 0.25
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = RKCMaterial(; model=NeoHooke())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.001 # shear parameter

    #===== homogeneous isotropic extension =====#
    F = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_ii_exact = 3 * params.K * (λ - 1)
    # test normal stresses
    for i in 1:3
        @test σ_pd_mean[i,i] ≈ σ_ii_exact rtol=1e-2
    end
    # test shear stresses
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    # should be nearly zero for isotropic extension, so we test that it is
    # at least 7 orders of magnitude smaller than normal stress
    @test all(σ_pd_vm_vec .< 1e-7 * σ_ii_exact)
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    @test all(isapprox.(σ_pd_hydro_vec, σ_ii_exact; rtol=1e-2))

    #===== pure shear =====#
    F = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_12_21_exact = params.μ * β
    # test affected shear components
    for (i, j) in ((1,2), (2,1))
        @test σ_pd_mean[i,j] ≈ σ_12_21_exact rtol=1e-3
    end
    # test normal stresses, should be at least small compared to shear stress
    for i in 1:3
        @test σ_pd_mean[i,i] < 0.005 * σ_12_21_exact
    end
    # test other shear stresses, should be zero
    for (i, j) in ((2,3), (3,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = params.μ * β * √3
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-3))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    # should be nearly zero for pure shear, therefore at least small
    # compared to shear stress
    @test all(isapprox.(σ_pd_hydro_vec, 0; atol=0.005σ_12_21_exact))

    #===== uniform deformation in x =====#
    F = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_11_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_22_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_33_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_22_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_11_exact - σ_22_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_11_exact + 2 * σ_22_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-2))

    #===== uniform deformation in y =====#
    F = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_22_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_33_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_11_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_22_exact - σ_11_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_22_exact + 2 * σ_11_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-2))

    #===== uniform deformation in z =====#
    F = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_33_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_22_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_33_exact - σ_11_22_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_33_exact + 2 * σ_11_22_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-2))
end

@testitem "Exported stress RKCMaterial NeoHookePenalty" setup=[StressExport] begin
    Δx = 0.25
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = RKCMaterial(; model=NeoHookePenalty())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.001 # shear parameter

    #===== homogeneous isotropic extension =====#
    F = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_ii_exact = 3 * params.K * (λ - 1)
    # test normal stresses
    for i in 1:3
        @test σ_pd_mean[i,i] ≈ σ_ii_exact rtol=1e-2
    end
    # test shear stresses
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    # should be nearly zero for isotropic extension, so we test that it is
    # at least 7 orders of magnitude smaller than normal stress
    @test all(σ_pd_vm_vec .< 1e-7 * σ_ii_exact)
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    @test all(isapprox.(σ_pd_hydro_vec, σ_ii_exact; rtol=1e-2))

    #===== pure shear =====#
    F = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_12_21_exact = params.μ * β
    # test affected shear components
    for (i, j) in ((1,2), (2,1))
        @test σ_pd_mean[i,j] ≈ σ_12_21_exact rtol=1e-3
    end
    # test normal stresses, should be at least small compared to shear stress
    for i in 1:3
        @test σ_pd_mean[i,i] < 0.005 * σ_12_21_exact
    end
    # test other shear stresses, should be zero
    for (i, j) in ((2,3), (3,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = params.μ * β * √3
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-3))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    # should be nearly zero for pure shear, therefore at least small
    # compared to shear stress
    @test all(isapprox.(σ_pd_hydro_vec, 0; atol=0.005σ_12_21_exact))

    #===== uniform deformation in x =====#
    F = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_11_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_22_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_33_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_22_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_11_exact - σ_22_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_11_exact + 2 * σ_22_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-2))

    #===== uniform deformation in y =====#
    F = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_22_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_33_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_11_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_22_exact - σ_11_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_22_exact + 2 * σ_11_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-2))

    #===== uniform deformation in z =====#
    F = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_33_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_22_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_33_exact - σ_11_22_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_33_exact + 2 * σ_11_22_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-2))
end


@testitem "Exported stress RKCRMaterial SaintVenantKirchhoff" setup=[StressExport] begin
    Δx = 0.25
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = RKCRMaterial(; model=SaintVenantKirchhoff())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.001 # shear parameter

    #===== homogeneous isotropic extension =====#
    F = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_ii_exact = 3 * params.K * (λ - 1)
    # test normal stresses
    for i in 1:3
        @test σ_pd_mean[i,i] ≈ σ_ii_exact rtol=1e-3
    end
    # test shear stresses
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    # should be nearly zero for isotropic extension, so we test that it is
    # at least 7 orders of magnitude smaller than normal stress
    @test all(σ_pd_vm_vec .< 1e-7 * σ_ii_exact)
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    @test all(isapprox.(σ_pd_hydro_vec, σ_ii_exact; rtol=1e-3))

    #===== pure shear =====#
    F = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_12_21_exact = params.μ * β
    # test affected shear components
    for (i, j) in ((1,2), (2,1))
        @test σ_pd_mean[i,j] ≈ σ_12_21_exact rtol=1e-3
    end
    # test normal stresses, should be at least small compared to shear stress
    for i in 1:3
        @test σ_pd_mean[i,i] < 0.005 * σ_12_21_exact
    end
    # test other shear stresses, should be zero
    for (i, j) in ((2,3), (3,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = params.μ * β * √3
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-3))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    # should be nearly zero for pure shear, therefore at least small
    # compared to shear stress
    @test all(isapprox.(σ_pd_hydro_vec, 0; atol=0.005σ_12_21_exact))

    #===== uniform deformation in x =====#
    F = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_11_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_22_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_33_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_22_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_11_exact - σ_22_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_11_exact + 2 * σ_22_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in y =====#
    F = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_22_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_33_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_11_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_22_exact - σ_11_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_22_exact + 2 * σ_11_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in z =====#
    F = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_33_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_22_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_33_exact - σ_11_22_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_33_exact + 2 * σ_11_22_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))
end

@testitem "Exported stress RKCRMaterial LinearElastic" setup=[StressExport] begin
    Δx = 0.25
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = RKCRMaterial(; model=LinearElastic())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.001 # shear parameter

    #===== homogeneous isotropic extension =====#
    F = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_ii_exact = 3 * params.K * (λ - 1)
    # test normal stresses
    for i in 1:3
        @test σ_pd_mean[i,i] ≈ σ_ii_exact rtol=1e-3
    end
    # test shear stresses
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    # should be nearly zero for isotropic extension, so we test that it is
    # at least 7 orders of magnitude smaller than normal stress
    @test all(σ_pd_vm_vec .< 1e-7 * σ_ii_exact)
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    @test all(isapprox.(σ_pd_hydro_vec, σ_ii_exact; rtol=1e-3))

    #===== pure shear =====#
    F = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_12_21_exact = params.μ * β
    # test affected shear components
    for (i, j) in ((1,2), (2,1))
        @test σ_pd_mean[i,j] ≈ σ_12_21_exact rtol=1e-3
    end
    # test normal stresses, should be at least small compared to shear stress
    for i in 1:3
        @test σ_pd_mean[i,i] < 0.005 * σ_12_21_exact
    end
    # test other shear stresses, should be zero
    for (i, j) in ((2,3), (3,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = params.μ * β * √3
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-3))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    # should be nearly zero for pure shear, therefore at least small
    # compared to shear stress
    @test all(isapprox.(σ_pd_hydro_vec, 0; atol=0.005σ_12_21_exact))

    #===== uniform deformation in x =====#
    F = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_11_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_22_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_33_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_22_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_11_exact - σ_22_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_11_exact + 2 * σ_22_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in y =====#
    F = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_22_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_33_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_33_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_11_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_22_exact - σ_11_33_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_22_exact + 2 * σ_11_33_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))

    #===== uniform deformation in z =====#
    F = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    #--- cauchy stress ---#
    σ_pd_vec = get_export_field(body, ts, F, :cauchy_stress)
    σ_pd_mean = mean_tensor(σ_pd_vec)
    σ_33_exact = (params.λ + 2 * params.μ) * (λ - 1)
    σ_11_22_exact = params.λ * (λ - 1)
    @test σ_pd_mean[1,1] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[2,2] ≈ σ_11_22_exact rtol=1e-2
    @test σ_pd_mean[3,3] ≈ σ_33_exact rtol=1e-2
    # test shear stresses should be nearly zero
    for (i, j) in ((1,2), (2,3), (3,1), (2,1), (3,2), (1,3))
        @test σ_pd_mean[i,j] ≈ 0.0 atol=1e-5
    end
    #--- von Mises stress ---#
    σ_pd_vm_vec = get_export_field(body, ts, F, :von_mises_stress)
    σ_vm_exact = abs(σ_33_exact - σ_11_22_exact)
    @test all(isapprox.(σ_pd_vm_vec, σ_vm_exact; rtol=1e-2))
    #--- hydrostatic stress ---#
    σ_pd_hydro_vec = get_export_field(body, ts, F, :hydrostatic_stress)
    σ_hydro_exact = (σ_33_exact + 2 * σ_11_22_exact) / 3
    @test all(isapprox.(σ_pd_hydro_vec, σ_hydro_exact; rtol=1e-3))
end
