@testitem "Internal force density correspondence" begin
    ref_position = [0.0 1.0 0.0 0.0 2.0
                    0.0 0.0 1.0 0.0 2.0
                    0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    δ = 1.5
    body = Body(NOSBMaterial(), ref_position, volume)
    material!(body, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0)
    failure_permit!(body, false)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; mat, storage, system, paramsetup) = chunk
    params = paramsetup
    (; position, b_int) = storage

    @test position == ref_position
    @test b_int == zeros(3, 5)

    # Boundary Condition:
    # Point 2 with v_z = 1 m/s with Δt = 0.0015 s
    position[1, 2] = 1.0015

    Peridynamics.calc_force_density!(chunk)

    @test b_int[:,1] ≈ [0.007185432432164514, 0.0023974081843325646, 0.0023974081843347313]
    @test b_int[:,2] ≈ [-0.007185432432123352, -9.047513544240966e-15, -6.0453612662025915e-15]
    @test b_int[:,3] ≈ [-2.077199268438824e-14, -0.002397408184361381, 3.2228109102914035e-14]
    @test b_int[:,4] ≈ [-2.0389697546832128e-14, 3.786380931014577e-14, -0.002397408184360914]
    @test b_int[:,5] ≈ [0.0, 0.0, 0.0]
end

@testitem "Internal force density correspondence interface" begin
    ref_position = [0.0 1.0 0.0 0.0 2.0
                    0.0 0.0 1.0 0.0 2.0
                    0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    δ = 1.5
    body = Body(NOSBMaterial(), ref_position, volume)
    point_set!(body, :a, [1])
    point_set!(body, :b, [2,3,4,5])
    material!(body, :a, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0)
    material!(body, :b, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0)
    failure_permit!(body, false)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; position, b_int) = chunk.storage

    @test position == ref_position
    @test b_int == zeros(3, 5)

    # Boundary Condition:
    # Point 2 with v_z = 1 m/s with Δt = 0.0015 s
    position[1, 2] = 1.0015

    Peridynamics.calc_force_density!(chunk)

    @test b_int[:,1] ≈ [0.007185432432164514, 0.0023974081843325646, 0.0023974081843347313]
    @test b_int[:,2] ≈ [-0.007185432432123352, -9.047513544240966e-15, -6.0453612662025915e-15]
    @test b_int[:,3] ≈ [-2.077199268438824e-14, -0.002397408184361381, 3.2228109102914035e-14]
    @test b_int[:,4] ≈ [-2.0389697546832128e-14, 3.786380931014577e-14, -0.002397408184360914]
    @test b_int[:,5] ≈ [0.0, 0.0, 0.0]
end
