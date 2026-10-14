# `CKIMaterial`: the continuum-kinematics-inspired formulation of
# `src/physics/continuum_kinematics_inspired.jl`.

@testitem "CKIMaterial: force density of a displaced point" begin
    ref_position = [0.0 1.0 0.0 0.0 2.0
                    0.0 0.0 1.0 0.0 2.0
                    0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    δ = 1.5
    body = Body(CKIMaterial(), ref_position, volume)
    @test_logs (:warn, r"specified manually") material!(body, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0, C1=1e11, C2=1e11, C3=1e11)
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

    @test b_int[:,1] ≈ [1.0000000000000625e9, 4.0060000000002503e8, 4.0060000000002503e8]
    @test b_int[:,2] ≈ [-1.449731150462047e9, 2.245287820579052e8, 2.245287820579052e8]
    @test b_int[:,3] ≈ [2.2486557523099208e8, -7.794353024117153e8, 1.543065203537848e8]
    @test b_int[:,4] ≈ [2.2486557523099208e8, 1.543065203537848e8, -7.794353024117153e8]
    @test b_int[:,5] ≈ [0.0, 0.0, 0.0]
end

@testitem "CKIMaterial: force density across a material interface" begin
    ref_position = [0.0 1.0 0.0 0.0 2.0
                    0.0 0.0 1.0 0.0 2.0
                    0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    δ = 1.5
    body = Body(CKIMaterial(), ref_position, volume)
    point_set!(body, :a, [1])
    point_set!(body, :b, [2,3,4,5])
    @test_logs (:warn, r"specified manually") material!(body, :a, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0, C1=1e11, C2=1e11, C3=1e11)
    @test_logs (:warn, r"specified manually") material!(body, :b, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0, C1=1e11, C2=1e11, C3=1e11)
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

    @test b_int[:,1] ≈ [1.0000000000000625e9, 4.0060000000002503e8, 4.0060000000002503e8]
    @test b_int[:,2] ≈ [-1.449731150462047e9, 2.245287820579052e8, 2.245287820579052e8]
    @test b_int[:,3] ≈ [2.2486557523099208e8, -7.794353024117153e8, 1.543065203537848e8]
    @test b_int[:,4] ≈ [2.2486557523099208e8, 1.543065203537848e8, -7.794353024117153e8]
    @test b_int[:,5] ≈ [0.0, 0.0, 0.0]
end

@testitem "standard_break_bond! / standard_break_bonds!: the InteractionSystem path" setup=[Fixtures] begin
    import Peridynamics: break_bond!, break_bonds!, each_one_ni_idx, bond_is_active

    body = Fixtures.cube(CKIMaterial(); n=4, m=2.015)
    chunk = Fixtures.chunk(body)
    (; storage, system, mat) = chunk
    dmg = Peridynamics.get_dmgmodel(mat)

    @test all(storage.one_ni_active)

    # breaking one one-neighbor interaction of point 1 flips only that flag
    one_ni_id = first(each_one_ni_idx(system, 1))
    break_bond!(storage, system, dmg, 1, one_ni_id)
    @test storage.one_ni_active[one_ni_id] == false
    @test bond_is_active(storage, system, one_ni_id) == false
    @test count(!, storage.one_ni_active) == 1

    # breaking all one-neighbor interactions of another point kills the whole point
    i = 2
    break_bonds!(storage, system, dmg, i)
    @test all(!bond_is_active(storage, system, id) for id in each_one_ni_idx(system, i))
    @test storage.n_active_one_nis[i] == 0
end
