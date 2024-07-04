@testitem "Energy surface correction bond based" begin
    ref_position = [0.0 1.0; 0.0 1.0; 0.0 1.0]
    volume = [1.0, 1.0]
    δ = 2.0
    E = 210e9
    body = Body(BBMaterial{EnergySurfaceCorrection}(), ref_position, volume)
    material!(body, horizon=δ, rho=1, E=E, Gc=1.0)
    failure_permit!(body, false)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; system) = chunk
    (; mfactor, scfactor) = system.correction

    @test mfactor == zeros(3, 2)
    @test scfactor == ones(2)

    Peridynamics.initialize!(dh, VelocityVerlet(steps=1))

    @test Peridynamics.analytical_stendens(E, 0.25, 0.001) ≈ 0.6 * 1e-6 * E
    @test mfactor ≈ fill(52.20262574164076, 3, 2)
    @test scfactor ≈ fill(52.20262574164076, 2)
end

@testitem "Energy surface correction bond based interface" begin
    ref_position = [0.0 1.0; 0.0 1.0; 0.0 1.0]
    volume = [1.0, 1.0]
    δ = 2.0
    Ea = 105e9
    Eb = 2 * Ea
    E = (Ea + Eb) / 2
    body = Body(BBMaterial{EnergySurfaceCorrection}(), ref_position, volume)
    point_set!(body, :a, [1])
    point_set!(body, :b, [2])
    material!(body, :a, horizon=δ, rho=1, E=Ea, Gc=1.0)
    material!(body, :b, horizon=δ, rho=1, E=Eb, Gc=1.0)
    failure_permit!(body, false)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; system) = chunk
    (; mfactor, scfactor) = system.correction

    @test mfactor == zeros(3, 2)
    @test scfactor == ones(2)

    Peridynamics.initialize!(dh, VelocityVerlet(steps=1))

    @test Peridynamics.analytical_stendens(E, 0.25, 0.001) ≈ 0.6 * 1e-6 * E
    @test mfactor ≈ fill(52.20262574164076, 3, 2)
    @test scfactor ≈ fill(52.20262574164076, 2)
end
