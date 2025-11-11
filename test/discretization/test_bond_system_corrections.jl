@testitem "EnergySurfaceCorrection BBMaterial" begin
    Δx = 0.1
    pos, vol = uniform_box(1, 1, 1, Δx)
    horizon = 3.01 * Δx
    rho = 8000
    E = 210e9
    nu = 0.25
    mat = BBMaterial{EnergySurfaceCorrection}()
    body = Body(mat, pos, vol)
    material!(body; horizon, rho, E, nu)
    ts = VelocityVerlet(steps=1)

    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; system) = chunk
    (; correction) = system
    (; mfactor, scfactor) = correction

    @test_broken minimum(mfactor[1, :]) ≈ 1.0 atol=0.08
    @test maximum(mfactor[1, :]) ≈ 3.7 atol=1.0
    @test_broken minimum(mfactor[2, :]) ≈ 1.0 atol=0.08
    @test maximum(mfactor[2, :]) ≈ 3.7 atol=1.0
    @test_broken minimum(mfactor[3, :]) ≈ 1.0 atol=0.08
    @test maximum(mfactor[3, :]) ≈ 3.7 atol=1.0
    @test_broken minimum(scfactor) ≈ 1.0 atol=0.08
    @test maximum(scfactor) ≈ 3.5 atol=1.0
end

@testitem "EnergySurfaceCorrection GBBMaterial" begin
    Δx = 0.1
    pos, vol = uniform_box(1, 1, 1, Δx)
    horizon = 3.01 * Δx
    rho = 8000
    E = 210e9
    nu = 0.25
    mat = GBBMaterial{EnergySurfaceCorrection}()
    body = Body(mat, pos, vol)
    material!(body; horizon, rho, E, nu)
    ts = VelocityVerlet(steps=1)

    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; system) = chunk
    (; correction) = system
    (; mfactor, scfactor) = correction

    @test minimum(mfactor[1, :]) ≈ 0.75 atol=0.2
    @test maximum(mfactor[1, :]) ≈ 1.5 atol=0.3
    @test minimum(mfactor[2, :]) ≈ 0.75 atol=0.2
    @test maximum(mfactor[2, :]) ≈ 1.5 atol=0.3
    @test minimum(mfactor[3, :]) ≈ 0.75 atol=0.2
    @test maximum(mfactor[3, :]) ≈ 1.5 atol=0.3
    @test minimum(scfactor) ≈ 0.75 atol=0.2
    @test maximum(scfactor) ≈ 1.5 atol=0.3
end

@testitem "EnergySurfaceCorrection OSBMaterial" begin
    Δx = 0.1
    pos, vol = uniform_box(1, 1, 1, Δx)
    horizon = 3.01 * Δx
    rho = 8000
    E = 210e9
    nu = 0.25
    mat = OSBMaterial{EnergySurfaceCorrection}()
    body = Body(mat, pos, vol)
    material!(body; horizon, rho, E, nu)
    ts = VelocityVerlet(steps=1)

    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; system) = chunk
    (; correction) = system
    (; mfactor, scfactor) = correction

    @test minimum(mfactor[1, :]) ≈ 0.75 atol=0.2
    @test maximum(mfactor[1, :]) ≈ 1.5 atol=0.3
    @test minimum(mfactor[2, :]) ≈ 0.75 atol=0.2
    @test maximum(mfactor[2, :]) ≈ 1.5 atol=0.3
    @test minimum(mfactor[3, :]) ≈ 0.75 atol=0.2
    @test maximum(mfactor[3, :]) ≈ 1.5 atol=0.3
    @test minimum(scfactor) ≈ 0.75 atol=0.2
    @test maximum(scfactor) ≈ 1.5 atol=0.3
end

@testitem "EnergySurfaceCorrection DHBBMaterial" begin
    Δx = 0.1
    pos, vol = uniform_box(1, 1, 1, Δx)
    horizon = 3.01 * Δx
    rho = 8000
    E = 210e9
    nu = 0.25
    mat = DHBBMaterial{EnergySurfaceCorrection}()
    body = Body(mat, pos, vol)
    material!(body; horizon, rho, E, nu)
    ts = VelocityVerlet(steps=1)

    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; system) = chunk
    (; correction) = system
    (; mfactor, scfactor) = correction

    @test_broken minimum(mfactor[1, :]) ≈ 1.0 atol=0.08
    @test maximum(mfactor[1, :]) ≈ 3.7 atol=1.0
    @test_broken minimum(mfactor[2, :]) ≈ 1.0 atol=0.08
    @test maximum(mfactor[2, :]) ≈ 3.7 atol=1.0
    @test_broken minimum(mfactor[3, :]) ≈ 1.0 atol=0.08
    @test maximum(mfactor[3, :]) ≈ 3.7 atol=1.0
    @test_broken minimum(scfactor) ≈ 1.0 atol=0.08
    @test maximum(scfactor) ≈ 3.5 atol=1.0
end

@testitem "Analytical strain energy density functions" begin
    test_cases = [
        (210e9, 0.25, 1.01),
        (200e9, 0.30, 1.05),
        (150e9, 0.20, 1.10),
        (300e9, 0.35, 1.02),
    ]
    for (E, nu, λ) in test_cases
        λ_lame = E * nu / ((1 + nu) * (1 - 2 * nu))
        μ_lame = E / (2 * (1 + nu))
        lame = [λ_lame, μ_lame]
        Ψ_small = 1/2 * λ_lame * (λ - 1)^2 + μ_lame * (λ - 1)^2
        Ψ_finite = 1/8 * λ_lame * (λ^2 - 1)^2 + 1/4 * μ_lame * (λ^2 - 1)^2
        @test Peridynamics.stendens_uniext_small_strain(lame, λ) ≈ Ψ_small
        @test Peridynamics.stendens_uniext_finite_strain(lame, λ) ≈ Ψ_finite
    end
end

@testitem "Get averaged Lamé parameters" begin
    pos = [0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
           0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    vol = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    horizon = 1.01
    rho = 8000
    E = 105e9
    nu = 0.25
    mat = BBMaterial{EnergySurfaceCorrection}()
    body = Body(mat, pos, vol)
    point_set!(x -> x < 4.5, body, :left)
    point_set!(x -> x > 4.5, body, :right)
    material!(body, :left; horizon, rho, E, nu)
    material!(body, :right; horizon, rho, E=2E, nu)
    ts = VelocityVerlet(steps=1)

    dh = Peridynamics.threads_data_handler(body, ts, 1)
    Peridynamics.initialize!(dh, ts)
    chunk = dh.chunks[1]
    (; system, storage, paramsetup) = chunk

    params1 = Peridynamics.get_params(paramsetup, 1)
    params4 = Peridynamics.get_params(paramsetup, 4)
    params5 = Peridynamics.get_params(paramsetup, 5)
    params6 = Peridynamics.get_params(paramsetup, 6)
    params7 = Peridynamics.get_params(paramsetup, 7)
    params10 = Peridynamics.get_params(paramsetup, 10)

    # point #1
    lame = Peridynamics.get_averaged_lame_parameters(system, storage, paramsetup, 1)
    @test lame[1] ≈ params1.λ
    @test lame[2] ≈ params1.μ

    # point #4
    lame = Peridynamics.get_averaged_lame_parameters(system, storage, paramsetup, 4)
    @test lame[1] ≈ params4.λ
    @test lame[2] ≈ params4.μ

    # point #5
    lame = Peridynamics.get_averaged_lame_parameters(system, storage, paramsetup, 5)
    @test lame[1] ≈ (params4.λ + 2 * params5.λ + params6.λ) / 4
    @test lame[2] ≈ (params4.μ + 2 * params5.μ + params6.μ) / 4

    # point #6
    lame = Peridynamics.get_averaged_lame_parameters(system, storage, paramsetup, 6)
    @test lame[1] ≈ (params5.λ + 2 * params6.λ + params7.λ) / 4
    @test lame[2] ≈ (params5.μ + 2 * params6.μ + params7.μ) / 4

    # point #7
    lame = Peridynamics.get_averaged_lame_parameters(system, storage, paramsetup, 7)
    @test lame[1] ≈ params7.λ
    @test lame[2] ≈ params7.μ

    # point #10
    lame = Peridynamics.get_averaged_lame_parameters(system, storage, paramsetup, 10)
    @test lame[1] ≈ params10.λ
    @test lame[2] ≈ params10.μ
end
