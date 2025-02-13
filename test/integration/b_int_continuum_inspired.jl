@testitem "Internal force density continuum inspired" begin
    ref_position = [0.0 1.0 0.0 0.0 2.0
                    0.0 0.0 1.0 0.0 2.0
                    0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    δ = 1.5
    body = Body(CKIMaterial(), ref_position, volume)
    material!(body, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0, C1=1e11, C2=1e11, C3=1e11)
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

    Peridynamics.calc_force_density!(chunk)

    @test b_int[:,1] ≈ [1.0000000000000625e9, 4.0060000000002503e8, 4.0060000000002503e8]
    @test b_int[:,2] ≈ [-1.449731150462047e9, 2.245287820579052e8, 2.245287820579052e8]
    @test b_int[:,3] ≈ [2.2486557523099208e8, -7.794353024117153e8, 1.543065203537848e8]
    @test b_int[:,4] ≈ [2.2486557523099208e8, 1.543065203537848e8, -7.794353024117153e8]
    @test b_int[:,5] ≈ [0.0, 0.0, 0.0]
end

@testitem "Internal force density continuum inspired interface" begin
    ref_position = [0.0 1.0 0.0 0.0 2.0
                    0.0 0.0 1.0 0.0 2.0
                    0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    δ = 1.5
    body = Body(CKIMaterial(), ref_position, volume)
    point_set!(body, :a, [1])
    point_set!(body, :b, [2,3,4,5])
    material!(body, :a, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0, C1=1e11, C2=1e11, C3=1e11)
    material!(body, :b, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0, C1=1e11, C2=1e11, C3=1e11)
    no_failure!(body)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; position, b_int) = chunk.storage

    @test position == ref_position
    @test b_int == zeros(3, 5)

    # Boundary Condition:
    # Point 2 with v_z = 1 m/s with Δt = 0.0015 s
    position[1, 2] = 1.0015

    Peridynamics.calc_force_density!(chunk)

    @test b_int[:,1] ≈ [1.0000000000000625e9, 4.0060000000002503e8, 4.0060000000002503e8]
    @test b_int[:,2] ≈ [-1.449731150462047e9, 2.245287820579052e8, 2.245287820579052e8]
    @test b_int[:,3] ≈ [2.2486557523099208e8, -7.794353024117153e8, 1.543065203537848e8]
    @test b_int[:,4] ≈ [2.2486557523099208e8, 1.543065203537848e8, -7.794353024117153e8]
    @test b_int[:,5] ≈ [0.0, 0.0, 0.0]
end
