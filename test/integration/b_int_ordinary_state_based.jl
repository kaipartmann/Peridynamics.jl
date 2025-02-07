@testitem "Internal force density ordinary state based" begin
    ref_position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    δ = 1.5
    E = 210e9
    body = Body(OSBMaterial(), ref_position, volume)
    material!(body, horizon=δ, rho=1, E=E, nu=0.25, Gc=1.0)
    failure_permit!(body, false)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; mat, storage, system, paramsetup) = chunk
    params = paramsetup
    (; position, b_int) = storage

    @test position == ref_position
    @test b_int == zeros(3, 2)

    # Boundary Condition:
    # Point 2 with v_z = 1 m/s with Δt = 0.0015 s
    position[1, 2] = 1.0015

    Peridynamics.calc_force_density!(chunk, 0, 0)

    b12 = 3.780000000000143e9
    @test b_int ≈ [b12 -b12; 0.0 0.0; 0.0 0.0]
end

@testitem "Internal force density ordinary state based interface" begin
    ref_position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 2.0]
    δ = 1.5
    Ea = 105e9
    Eb = 2 * Ea
    E = (Ea + Eb) / 2
    body = Body(OSBMaterial(), ref_position, volume)
    point_set!(body, :a, [1])
    point_set!(body, :b, [2])
    material!(body, :a, horizon=δ, rho=1, E=Ea, nu=0.25, Gc=1.0)
    material!(body, :b, horizon=δ, rho=1, E=Eb, nu=0.25, Gc=1.0)
    failure_permit!(body, false)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; position, b_int) = chunk.storage

    @test position == ref_position
    @test b_int == zeros(3, 2)

    # Boundary Condition:
    # Point 2 with v_z = 1 m/s with Δt = 0.0015 s
    position[1, 2] = 1.0015

    Peridynamics.calc_force_density!(chunk, 0, 0)

    b1 = 4.252500000000161e9
    b2 = 2.1262500000000806e9
    @test b_int ≈ [b1 -b2; 0.0 0.0; 0.0 0.0]
end
