@testitem "Internal force density correspondence ZEMSilling" begin
    ref_position = [0.0 1.0 0.0 0.0 2.0
                    0.0 0.0 1.0 0.0 2.0
                    0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    δ = 1.5
    body = Body(CMaterial(model=MooneyRivlin()), ref_position, volume)
    material!(body, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0)
    no_failure!(body)

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

    Peridynamics.calc_force_density!(chunk, 0, 0)

    @test b_int[:,1] ≈ [0.007185432432164514, 0.0023974081843325646, 0.0023974081843347313]
    @test b_int[:,2] ≈ [-0.007185432432123352, -9.047513544240966e-15, -6.0453612662025915e-15]
    @test b_int[:,3] ≈ [-2.077199268438824e-14, -0.002397408184361381, 3.2228109102914035e-14]
    @test b_int[:,4] ≈ [-2.0389697546832128e-14, 3.786380931014577e-14, -0.002397408184360914]
    @test b_int[:,5] ≈ [0.0, 0.0, 0.0]
end

@testitem "Internal force density correspondence ZEMWan" begin
    ref_position = [0.0 1.0 0.0 0.0 2.0
                    0.0 0.0 1.0 0.0 2.0
                    0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    δ = 1.5
    body = Body(CMaterial(model=MooneyRivlin(), zem=ZEMWan()), ref_position, volume)
    material!(body, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0)
    no_failure!(body)

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

    Peridynamics.calc_force_density!(chunk, 0, 0)

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
    body = Body(CMaterial(model=MooneyRivlin()), ref_position, volume)
    point_set!(body, :a, [1])
    point_set!(body, :b, [2,3,4,5])
    material!(body, :a, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0)
    material!(body, :b, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0)
    no_failure!(body)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; position, b_int) = chunk.storage

    @test position == ref_position
    @test b_int == zeros(3, 5)

    # Boundary Condition:
    # Point 2 with v_z = 1 m/s with Δt = 0.0015 s
    position[1, 2] = 1.0015

    Peridynamics.calc_force_density!(chunk, 0, 0)

    @test b_int[:,1] ≈ [0.007185432432164514, 0.0023974081843325646, 0.0023974081843347313]
    @test b_int[:,2] ≈ [-0.007185432432123352, -9.047513544240966e-15, -6.0453612662025915e-15]
    @test b_int[:,3] ≈ [-2.077199268438824e-14, -0.002397408184361381, 3.2228109102914035e-14]
    @test b_int[:,4] ≈ [-2.0389697546832128e-14, 3.786380931014577e-14, -0.002397408184360914]
    @test b_int[:,5] ≈ [0.0, 0.0, 0.0]
end

@testitem "Internal force density correspondence rotated" begin
    ref_position = [0.0 1.0 0.0 0.0 2.0
                    0.0 0.0 1.0 0.0 2.0
                    0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    δ = 1.5
    body = Body(CRMaterial(model=LinearElastic()), ref_position, volume)
    material!(body, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0)
    no_failure!(body)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; mat, storage, system, paramsetup) = chunk
    params = paramsetup
    (; position, b_int) = storage

    @test position == ref_position
    @test b_int == zeros(3, 5)

    # Boundary Condition:
    position[1, 2] = 1.1

    Peridynamics.calc_force_density!(chunk, 1e-7, 1e-7)

    @test isapprox(b_int, [
        -3.85668e-16  -2.9442e-14    1.49138e-14   1.49138e-14  0.0
         1.34581e-14   9.0082e-15   -2.36931e-14   1.22675e-15  0.0
         1.59116e-14   1.39152e-14  -1.22675e-15  -2.86001e-14  0.0
    ]; atol=3eps())
end

@testitem "Internal force density correspondence rotated ZEMWan" begin
    ref_position = [0.0 1.0 0.0 0.0 2.0
                    0.0 0.0 1.0 0.0 2.0
                    0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    δ = 1.5
    body = Body(CRMaterial(model=LinearElastic(), zem=ZEMWan()), ref_position, volume)
    material!(body, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0)
    no_failure!(body)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; mat, storage, system, paramsetup) = chunk
    params = paramsetup
    (; position, b_int) = storage

    @test position == ref_position
    @test b_int == zeros(3, 5)

    # Boundary Condition:
    position[1, 2] = 1.1

    Peridynamics.calc_force_density!(chunk, 1e-7, 1e-7)

    @test isapprox(b_int, [
        6.99192e-16  -2.71841e-15   1.04101e-15   9.78208e-16  0.0
        1.73062e-15   6.53589e-16  -2.2012e-15   -1.83013e-16  0.0
        2.00246e-15   1.01964e-15  -3.9205e-16   -2.63005e-15  0.0
    ]; atol=3eps())
end
