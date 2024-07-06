@testitem "Internal force density bond based" begin
    ref_position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    δ = 1.5
    E = 210e9
    body = Body(BBMaterial(), ref_position, volume)
    material!(body, horizon=δ, rho=1, E=E, Gc=1.0)
    failure_permit!(body, false)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; position, b_int) = chunk.storage
    (; bonds) = chunk.system

    @test position == ref_position
    @test b_int == zeros(3, 2)

    # Boundary Condition:
    # Point 2 with v_z = 1 m/s with Δt = 0.0015 s
    position[1, 2] = 1.0015

    Peridynamics.calc_force_density!(chunk)

    b12 = 18 * E / (3 * (1 - 2 * 0.25)) / (π * δ^4) * 1.0015 * 0.0015/1.0015 * 1.0
    @test b_int ≈ [b12 -b12; 0.0 0.0; 0.0 0.0]
end

@testitem "Internal force density bond based interface" begin
    ref_position = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    volume = [1.0, 1.0]
    δ = 1.5
    Ea = 105e9
    Eb = 2 * Ea
    E = (Ea + Eb) / 2
    body = Body(BBMaterial(), ref_position, volume)
    point_set!(body, :a, [1])
    point_set!(body, :b, [2])
    material!(body, :a, horizon=δ, rho=1, E=Ea, Gc=1.0)
    material!(body, :b, horizon=δ, rho=1, E=Eb, Gc=1.0)
    failure_permit!(body, false)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; position, b_int) = chunk.storage
    (; bonds) = chunk.system

    @test position == ref_position
    @test b_int == zeros(3, 2)

    # Boundary Condition:
    # Point 2 with v_z = 1 m/s with Δt = 0.0015 s
    position[1, 2] = 1.0015

    Peridynamics.calc_force_density!(chunk)

    b12 = 18 * E / (3 * (1 - 2 * 0.25)) / (π * δ^4) * 1.0015 * 0.0015/1.0015 * 1.0
    @test b_int ≈ [b12 -b12; 0.0 0.0; 0.0 0.0]
end
