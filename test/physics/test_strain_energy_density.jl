@testsnippet PsiExport begin
    using Peridynamics.StaticArrays, Peridynamics.Printf
    mean(x) = sum(x) / length(x)
    function test_stendens(body, ts, F_a, Ψ_a, tols; testcase="", atol0=eps())
        dh = Peridynamics.threads_data_handler(body, ts, 1)
        Peridynamics.initialize!(dh, ts)
        (; mat, system, storage, paramsetup) = dh.chunks[1]
        field = Val{:strain_energy_density}()

        # no deformation -> no strain energy density
        storage.position .= system.position
        storage.strain_energy_density .= 0.0
        Ψ_0 = Peridynamics.export_field(field, mat, system, storage, paramsetup, 0.0)
        @test all(isapprox.(Ψ_0, 0.0; atol=atol0))

        # apply deformation gradient F_a
        for i in Peridynamics.each_point_idx(system)
            Xi = Peridynamics.get_vector(system.position, i)
            xi = F_a * Xi
            Peridynamics.update_vector!(storage.position, i, xi)
        end
        Peridynamics.calc_force_density!(dh, 0.0, 0.0)

        Ψ_pd = Peridynamics.export_field(field, mat, system, storage, paramsetup, 0.0)
        Ψ̂_pd = mean(Ψ_pd)
        ΔΨ = (Ψ_pd .- Ψ_a) ./ Ψ_a
        ΔΨ_max, ΔΨ_min = maximum(ΔΨ), minimum(ΔΨ)

        # write detailed printed output that shows the calculation results
        println("-"^20 * "  Ψ TEST  " * "-"^20)
        if !isempty(testcase)
            @printf("Test case: %s\n", testcase)
        end
        @printf("material: %s\n", nameof(typeof(mat)))
        if hasproperty(mat, :constitutive_model)
            @printf("model: %s\n", mat.constitutive_model)
        end
        @printf("Expected Ψ:        %.6g\n", Ψ_a)
        @printf("Calc. mean(Ψ):     %.6g\n", Ψ̂_pd)
        @printf("  mean error:      %7.2f %%\n", (Ψ̂_pd - Ψ_a) / Ψ_a * 100)
        @printf("  min rel. error:  %7.2f %%\n", ΔΨ_min * 100)
        @printf("  max rel. error:  %7.2f %%\n", ΔΨ_max * 100)

        @test ΔΨ_min > tols[1]
        @test ΔΨ_max < tols[2]
        return nothing
    end
end

@testitem "Strain energy density export BBMaterial" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    pos, vol = uniform_box(1,1,1,Δx)
    mat = BBMaterial{NoCorrection}()
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    testcase = "homogeneous isotropic extension"
    λ = 1.1
    ε = λ - 1
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/2 * params.K * ε^2
    tols = (-0.9, 0.3)
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "pure shear deformation"
    β = 0.1
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = 1/2 * params.G * β^2
    tols = (-0.9, 0.3)
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in x-direction"
    λ = 1.1
    ε = λ - 1
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = 3/5 * params.E * ε^2
    tols = (-0.9, 0.3)
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in y-direction"
    λ = 1.1
    ε = λ - 1
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    Ψ_a = 3/5 * params.E * ε^2
    tols = (-0.9, 0.3)
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in z-direction"
    λ = 1.1
    ε = λ - 1
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    Ψ_a = 3/5 * params.E * ε^2
    tols = (-0.9, 0.3)
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)
end

@testitem "Strain energy density export DHBBMaterial" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    pos, vol = uniform_box(1,1,1,Δx)
    mat = DHBBMaterial{NoCorrection}()
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    testcase = "homogeneous isotropic extension"
    λ = 1.1
    ε = λ - 1
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/2 * params.K * ε^2
    tols = (-0.9, 0.3)
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "pure shear deformation"
    β = 0.1
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = 1/2 * params.G * β^2
    tols = (-0.9, 0.3)
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in x-direction"
    λ = 1.1
    ε = λ - 1
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = 3/5 * params.E * ε^2
    tols = (-0.9, 0.3)
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in y-direction"
    λ = 1.1
    ε = λ - 1
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    Ψ_a = 3/5 * params.E * ε^2
    tols = (-0.9, 0.3)
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in z-direction"
    λ = 1.1
    ε = λ - 1
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    Ψ_a = 3/5 * params.E * ε^2
    tols = (-0.9, 0.3)
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)
end

@testitem "Strain energy density export OSBMaterial" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = OSBMaterial{NoCorrection}()
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    testcase = "homogeneous isotropic extension"
    λ = 1.1
    ε = λ - 1
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/2 * params.K * ε^2
    tols = (-1e-10, 1e-10) # should be very accurate for a homogeneous iso. extension
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "pure shear deformation"
    β = 0.1
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = 1/2 * params.G * β^2
    tols = (-0.4, 0.4) # higher errors
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in x-direction"
    λ = 1.1
    ε = λ - 1
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = 1/2 * params.λ * ε^2 + params.μ * ε^2
    tols = (-0.9, 0.3)
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in y-direction"
    λ = 1.1
    ε = λ - 1
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    Ψ_a = 1/2 * params.λ * ε^2 + params.μ * ε^2
    tols = (-0.9, 0.3)
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in z-direction"
    λ = 1.1
    ε = λ - 1
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    Ψ_a = 1/2 * params.λ * ε^2 + params.μ * ε^2
    tols = (-0.9, 0.3)
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)
end

@testitem "Strain energy density export CMaterial LinearElastic" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = CMaterial(; model=LinearElastic())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    # all tests should have very high accuracy
    tols = (-1e-12, 1e-12)
    λ = 1.1 # uniform extension
    ε = λ - 1 # strain
    β = 0.1 # shear parameter

    testcase = "homogeneous isotropic extension"
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/8 * params.K * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "pure shear deformation"
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = (params.μ/2) * β^2 + (params.λ/8 + params.μ/4) * β^4
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in x-direction"
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = (params.λ + 2 * params.μ)/8 * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in y-direction"
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction

    testcase = "uniform extension in z-direction"
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction
end

@testitem "Strain energy density export CMaterial SaintVenantKirchhoff" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = CMaterial(; model=SaintVenantKirchhoff())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    # all tests should have very high accuracy
    tols = (-1e-12, 1e-12)
    λ = 1.1 # uniform extension
    ε = λ - 1 # strain
    β = 0.1 # shear parameter

    testcase = "homogeneous isotropic extension"
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/8 * params.K * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "pure shear deformation"
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = (params.μ/2) * β^2 + (params.λ/8 + params.μ/4) * β^4
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in x-direction"
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = (params.λ + 2 * params.μ)/8 * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in y-direction"
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction

    testcase = "uniform extension in z-direction"
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction
end

@testitem "Strain energy density export CMaterial NeoHooke" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = CMaterial(; model=NeoHooke())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    tols = (-0.02, 0.02) # slightly less accurate due to nonlinearity
    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.05 # shear parameter

    testcase = "homogeneous isotropic extension"
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/8 * params.K * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4)

    testcase = "pure shear deformation"
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = (params.μ/2) * β^2 + (params.λ/8 + params.μ/4) * β^4
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4)

    testcase = "uniform extension in x-direction"
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = (params.λ + 2 * params.μ)/8 * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4)

    testcase = "uniform extension in y-direction"
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4) # same Ψ_a as x-direction

    testcase = "uniform extension in z-direction"
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4) # same Ψ_a as x-direction
end

@testitem "Strain energy density export CMaterial NeoHookePenalty" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = CMaterial(; model=NeoHookePenalty())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    tols = (-0.02, 0.02) # slightly less accurate due to nonlinearity
    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.05 # shear parameter

    testcase = "homogeneous isotropic extension"
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/8 * params.K * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4)

    testcase = "pure shear deformation"
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = (params.μ/2) * β^2 + (params.λ/8 + params.μ/4) * β^4
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4)

    testcase = "uniform extension in x-direction"
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = (params.λ + 2 * params.μ)/8 * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4)

    testcase = "uniform extension in y-direction"
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4) # same Ψ_a as x-direction

    testcase = "uniform extension in z-direction"
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4) # same Ψ_a as x-direction
end

@testitem "Strain energy density export CRMaterial LinearElastic" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = CRMaterial(; model=LinearElastic())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    # all tests should have very high accuracy
    tols = (-1e-12, 1e-12)
    λ = 1.1 # uniform extension
    ε = λ - 1 # strain
    β = 0.1 # shear parameter

    testcase = "homogeneous isotropic extension"
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/8 * params.K * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "pure shear deformation"
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = (params.μ/2) * β^2 + (params.λ/8 + params.μ/4) * β^4
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in x-direction"
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = (params.λ + 2 * params.μ)/8 * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in y-direction"
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction

    testcase = "uniform extension in z-direction"
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction
end

@testitem "Strain energy density export CRMaterial SaintVenantKirchhoff" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = CRMaterial(; model=SaintVenantKirchhoff())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    # all tests should have very high accuracy
    tols = (-1e-12, 1e-12)
    λ = 1.1 # uniform extension
    ε = λ - 1 # strain
    β = 0.1 # shear parameter

    testcase = "homogeneous isotropic extension"
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/8 * params.K * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "pure shear deformation"
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = (params.μ/2) * β^2 + (params.λ/8 + params.μ/4) * β^4
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in x-direction"
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = (params.λ + 2 * params.μ)/8 * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in y-direction"
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction

    testcase = "uniform extension in z-direction"
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction
end

@testitem "Strain energy density export RKCMaterial LinearElastic" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = RKCMaterial(; model=LinearElastic())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    # all tests should have very high accuracy
    tols = (-1e-12, 1e-12)
    λ = 1.1 # uniform extension
    ε = λ - 1 # strain
    β = 0.1 # shear parameter

    testcase = "homogeneous isotropic extension"
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/8 * params.K * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "pure shear deformation"
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = (params.μ/2) * β^2 + (params.λ/8 + params.μ/4) * β^4
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in x-direction"
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = (params.λ + 2 * params.μ)/8 * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in y-direction"
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction

    testcase = "uniform extension in z-direction"
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction
end

@testitem "Strain energy density export RKCMaterial SaintVenantKirchhoff" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = RKCMaterial(; model=SaintVenantKirchhoff())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    # all tests should have very high accuracy
    tols = (-1e-12, 1e-12)
    λ = 1.1 # uniform extension
    ε = λ - 1 # strain
    β = 0.1 # shear parameter

    testcase = "homogeneous isotropic extension"
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/8 * params.K * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "pure shear deformation"
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = (params.μ/2) * β^2 + (params.λ/8 + params.μ/4) * β^4
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in x-direction"
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = (params.λ + 2 * params.μ)/8 * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in y-direction"
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction

    testcase = "uniform extension in z-direction"
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction
end

@testitem "Strain energy density export RKCMaterial NeoHooke" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = RKCMaterial(; model=NeoHooke())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    tols = (-0.02, 0.02) # slightly less accurate due to nonlinearity
    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.05 # shear parameter

    testcase = "homogeneous isotropic extension"
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/8 * params.K * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4)

    testcase = "pure shear deformation"
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = (params.μ/2) * β^2 + (params.λ/8 + params.μ/4) * β^4
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4)

    testcase = "uniform extension in x-direction"
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = (params.λ + 2 * params.μ)/8 * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4)

    testcase = "uniform extension in y-direction"
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4) # same Ψ_a as x-direction

    testcase = "uniform extension in z-direction"
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4) # same Ψ_a as x-direction
end

@testitem "Strain energy density export RKCMaterial NeoHookePenalty" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = RKCMaterial(; model=NeoHookePenalty())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    tols = (-0.02, 0.02) # slightly less accurate due to nonlinearity
    λ = 1.001 # uniform extension
    ε = λ - 1 # strain
    β = 0.05 # shear parameter

    testcase = "homogeneous isotropic extension"
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/8 * params.K * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4)

    testcase = "pure shear deformation"
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = (params.μ/2) * β^2 + (params.λ/8 + params.μ/4) * β^4
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4)

    testcase = "uniform extension in x-direction"
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = (params.λ + 2 * params.μ)/8 * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4)

    testcase = "uniform extension in y-direction"
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4) # same Ψ_a as x-direction

    testcase = "uniform extension in z-direction"
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase, atol0=1e-4) # same Ψ_a as x-direction
end

@testitem "Strain energy density export RKCRMaterial LinearElastic" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = RKCRMaterial(; model=LinearElastic())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    # all tests should have very high accuracy
    tols = (-1e-12, 1e-12)
    λ = 1.1 # uniform extension
    ε = λ - 1 # strain
    β = 0.1 # shear parameter

    testcase = "homogeneous isotropic extension"
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/8 * params.K * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "pure shear deformation"
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = (params.μ/2) * β^2 + (params.λ/8 + params.μ/4) * β^4
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in x-direction"
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = (params.λ + 2 * params.μ)/8 * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in y-direction"
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction

    testcase = "uniform extension in z-direction"
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction
end

@testitem "Strain energy density export RKCRMaterial SaintVenantKirchhoff" setup=[PsiExport] begin
    Δx = 0.2
    horizon = 3.01Δx
    E = 210e9
    nu = 0.25
    pos, vol = uniform_box(1,1,1,Δx)
    mat = RKCRMaterial(; model=SaintVenantKirchhoff())
    body = Body(mat, pos, vol)
    material!(body; horizon, rho=8000, E, nu)
    params = body.point_params[1]
    ts = VelocityVerlet(steps=1)

    # all tests should have very high accuracy
    tols = (-1e-12, 1e-12)
    λ = 1.1 # uniform extension
    ε = λ - 1 # strain
    β = 0.1 # shear parameter

    testcase = "homogeneous isotropic extension"
    F_a = @SMatrix [λ 0 0; 0 λ 0; 0 0 λ]
    Ψ_a = 9/8 * params.K * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "pure shear deformation"
    F_a = @SMatrix [1 β 0; 0 1 0; 0 0 1]
    Ψ_a = (params.μ/2) * β^2 + (params.λ/8 + params.μ/4) * β^4
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in x-direction"
    F_a = @SMatrix [λ 0 0; 0 1 0; 0 0 1]
    Ψ_a = (params.λ + 2 * params.μ)/8 * (λ^2 - 1)^2
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase)

    testcase = "uniform extension in y-direction"
    F_a = @SMatrix [1 0 0; 0 λ 0; 0 0 1]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction

    testcase = "uniform extension in z-direction"
    F_a = @SMatrix [1 0 0; 0 1 0; 0 0 λ]
    test_stendens(body, ts, F_a, Ψ_a, tols; testcase) # same Ψ_a as x-direction
end
