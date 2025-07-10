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
    # Point 2 with v_z = 1 m/s with Δt = 0.0015 s
    position[1, 2] = 1.0015

    Peridynamics.calc_force_density!(chunk, 0, 0)

    @test b_int[:,1] ≈ [0.0, 1.1846541080217141e-14, 1.59116313101076e-14]
    @test b_int[:,2] ≈ [-9.81399596527191e-15, 1.0619791584558151e-14, 1.3915198328621624e-14]
    @test b_int[:,3] ≈ [4.8349988789970084e-15, -2.3693082160434282e-14, -1.2267494956589887e-15]
    @test b_int[:,4] ≈ [4.978997086274901e-15, 1.2267494956589887e-15, -2.8600080143070235e-14]
    @test b_int[:,5] ≈ [0.0, 0.0, 0.0]
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
    # Point 2 with v_z = 1 m/s with Δt = 0.0015 s
    position[1, 2] = 1.0015

    Peridynamics.calc_force_density!(chunk, 0, 0)

    @test b_int[:,1] ≈ [7.105427357601002e-16, 1.5206787245922409e-15, 1.9710800241257256e-15]
    @test b_int[:,2] ≈ [-5.329070518200751e-16, 1.3758314660245038e-15, 1.5633181769341484e-15]
    @test b_int[:,3] ≈ [2.1581533343596185e-17, -2.4591894988755905e-15, -6.50044234103052e-16]
    @test b_int[:,4] ≈ [-1.992172172836213e-16, -4.373206917411543e-16, -2.884353966956822e-15]
    @test b_int[:,5] ≈ [0.0, 0.0, 0.0]
end
