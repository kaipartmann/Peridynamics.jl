@testitem "one- and two-neighbor interactions" begin
    position = [0.0 1.0 0.0
                0.0 0.0 0.0
                0.0 0.0 1.0]
    volume = fill(1.0, 3)
    body = Body(CKIMaterial(), position, volume)
    @test_logs (:warn, r"specified manually") material!(body; horizon=1.5, rho=8e-6, E=210e3, nu=0.3, Gc=1.0, C1=1, C2=1)
    pd = Peridynamics.PointDecomposition(body, 1)

    system = Peridynamics.InteractionSystem(body, pd, 1)

    @test Peridynamics.has_three_nis(body.point_params[1]) == false
    @test Peridynamics.has_three_nis(body,1) == false

    @test Peridynamics.has_two_nis(body,1) == true

    @test system.one_nis == [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(3, √2, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(2, √2, true)]
    @test system.n_one_nis == [2, 2, 2]
    @test system.one_ni_idxs == [1:2, 3:4, 5:6]
    @test system.two_nis == [
        Peridynamics.TwoNeighborInteraction(1, 2, 1.0),
        Peridynamics.TwoNeighborInteraction(3, 4, 1.0),
        Peridynamics.TwoNeighborInteraction(5, 6, 1.0)]
    @test system.n_two_nis == [1, 1, 1]
    @test system.two_ni_idxs == [1:1, 2:2, 3:3]
    @test system.three_nis == Vector{Peridynamics.ThreeNeighborInteraction}()
    @test system.n_three_nis == Vector{Int}()
    @test system.three_ni_idxs == Vector{UnitRange{Int}}()

    point_set!(body, :a, 1:2)
    @test_logs (:warn, r"specified manually") material!(body, :a; horizon=1.5, rho=7e-6, E=200e3, nu=0.3, Gc=1.0, C1=1, C2=1)

    ts = VelocityVerlet(steps=1)
    dh = Peridynamics.threads_data_handler(body, ts, 1)
    chunk = dh.chunks[1]
    @test Peridynamics.has_three_nis(chunk.paramsetup) == false
end

@testitem "one- and three-neighbor interactions" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 0.0 1.0
                0.0 0.0 1.0 0.0]
    volume = fill(1.0, 4)
    body = Body(CKIMaterial(), position, volume)
    @test_logs (:warn, r"specified manually") material!(body; horizon=1.5, rho=8e-6, E=210e3, nu=0.3, Gc=1.0, C1=1, C3=1)
    pd = Peridynamics.PointDecomposition(body, 1)

    system = Peridynamics.InteractionSystem(body, pd, 1)

    @test Peridynamics.has_two_nis(body.point_params[1]) == false
    @test Peridynamics.has_two_nis(body) == false
    @test Peridynamics.has_two_nis(body,1) == false

    @test Peridynamics.has_three_nis(body,1) == true

    @test system.one_nis == [
    Peridynamics.Bond(2, 1.0, true),
    Peridynamics.Bond(3, 1.0, true),
    Peridynamics.Bond(4, 1.0, true),
    Peridynamics.Bond(1, 1.0, true),
    Peridynamics.Bond(3, √2, true),
    Peridynamics.Bond(4, √2, true),
    Peridynamics.Bond(1, 1.0, true),
    Peridynamics.Bond(2, √2, true),
    Peridynamics.Bond(4, √2, true),
    Peridynamics.Bond(1, 1.0, true),
    Peridynamics.Bond(2, √2, true),
    Peridynamics.Bond(3, √2, true)]
    @test system.n_one_nis == [3, 3, 3, 3]
    @test system.one_ni_idxs == [1:3, 4:6, 7:9, 10:12]
    @test system.two_nis == Vector{Peridynamics.TwoNeighborInteraction}()
    @test system.n_two_nis == Vector{Int}()
    @test system.two_ni_idxs == Vector{UnitRange{Int}}()
    @test system.three_nis == [
        Peridynamics.ThreeNeighborInteraction(1, 2, 3, 1.0),
        Peridynamics.ThreeNeighborInteraction(4, 5, 6, 1.0),
        Peridynamics.ThreeNeighborInteraction(7, 8, 9, 1.0),
        Peridynamics.ThreeNeighborInteraction(10, 11, 12, 1.0)]
    @test system.n_three_nis == [1, 1, 1, 1]
    @test system.three_ni_idxs == [1:1, 2:2, 3:3, 4:4]

    point_set!(body, :a, 1:2)
    @test_logs (:warn, r"specified manually") material!(body, :a; horizon=1.5, rho=7e-6, E=200e3, nu=0.3, Gc=1.0, C1=1, C3=1)

    ts = VelocityVerlet(steps=1)
    dh = Peridynamics.threads_data_handler(body, ts, 1)
    chunk = dh.chunks[1]
    @test Peridynamics.has_two_nis(chunk.paramsetup) == false
end

@testitem "one-, two-, and three-neighbor interactions" begin
    position = [0.0 1.0 0.0 0.0 2.0
                0.0 0.0 1.0 0.0 2.0
                0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    body = Body(CKIMaterial(), position, volume)
    @test_logs (:warn, r"specified manually") material!(body; horizon=1.5, rho=8e-6, E=210e3, nu=0.3, Gc=1.0, C1=1, C2=1, C3=1)
    pd = Peridynamics.PointDecomposition(body, 1)

    system = Peridynamics.InteractionSystem(body, pd, 1)

    @test system.one_nis == [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(3, √2, true),
        Peridynamics.Bond(4, √2, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(2, √2, true),
        Peridynamics.Bond(4, √2, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(2, √2, true),
        Peridynamics.Bond(3, √2, true)]
    @test system.n_one_nis == [3, 3, 3, 3, 0]
    @test system.one_ni_idxs == [1:3, 4:6, 7:9, 10:12, 13:12]
    @test system.two_nis == [
        Peridynamics.TwoNeighborInteraction(1, 2, 1.0),
        Peridynamics.TwoNeighborInteraction(1, 3, 1.0),
        Peridynamics.TwoNeighborInteraction(2, 3, 1.0),
        Peridynamics.TwoNeighborInteraction(4, 5, 1.0),
        Peridynamics.TwoNeighborInteraction(4, 6, 1.0),
        Peridynamics.TwoNeighborInteraction(5, 6, √3),
        Peridynamics.TwoNeighborInteraction(7, 8, 1.0),
        Peridynamics.TwoNeighborInteraction(7, 9, 1.0),
        Peridynamics.TwoNeighborInteraction(8, 9, √3),
        Peridynamics.TwoNeighborInteraction(10, 11, 1.0),
        Peridynamics.TwoNeighborInteraction(10, 12, 1.0),
        Peridynamics.TwoNeighborInteraction(11, 12, √3)]
    @test system.n_two_nis == [3, 3, 3, 3, 0]
    @test system.two_ni_idxs == [1:3, 4:6, 7:9, 10:12, 0:-1]
    @test system.three_nis == [
        Peridynamics.ThreeNeighborInteraction(1, 2, 3, 1.0),
        Peridynamics.ThreeNeighborInteraction(4, 5, 6, 1.0),
        Peridynamics.ThreeNeighborInteraction(7, 8, 9, 1.0),
        Peridynamics.ThreeNeighborInteraction(10, 11, 12, 1.0)]
    @test system.n_three_nis == [1, 1, 1, 1, 0]
    @test system.three_ni_idxs == [1:1, 2:2, 3:3, 4:4, 0:-1]
end

@testitem "initialized InteractionSystem" begin
    position = [0.0 1.0 0.0 0.0 2.0
                0.0 0.0 1.0 0.0 2.0
                0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    body = Body(CKIMaterial(), position, volume)
    @test_logs (:warn, r"specified manually") material!(body; horizon=1.5, rho=8e-6, E=210e3, nu=0.3, Gc=1.0, C1=1, C2=1, C3=1)

    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    is = dh.chunks[1].system

    @test is.one_nis == [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(4, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(3, √2, true),
        Peridynamics.Bond(4, √2, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(2, √2, true),
        Peridynamics.Bond(4, √2, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(2, √2, true),
        Peridynamics.Bond(3, √2, true)]
    @test is.n_one_nis == [3, 3, 3, 3, 0]
    @test is.one_ni_idxs == [1:3, 4:6, 7:9, 10:12, 13:12]
    @test is.two_nis == [
        Peridynamics.TwoNeighborInteraction(1, 2, 1.0),
        Peridynamics.TwoNeighborInteraction(1, 3, 1.0),
        Peridynamics.TwoNeighborInteraction(2, 3, 1.0),
        Peridynamics.TwoNeighborInteraction(4, 5, 1.0),
        Peridynamics.TwoNeighborInteraction(4, 6, 1.0),
        Peridynamics.TwoNeighborInteraction(5, 6, √3),
        Peridynamics.TwoNeighborInteraction(7, 8, 1.0),
        Peridynamics.TwoNeighborInteraction(7, 9, 1.0),
        Peridynamics.TwoNeighborInteraction(8, 9, √3),
        Peridynamics.TwoNeighborInteraction(10, 11, 1.0),
        Peridynamics.TwoNeighborInteraction(10, 12, 1.0),
        Peridynamics.TwoNeighborInteraction(11, 12, √3)]
    @test is.n_two_nis == [3, 3, 3, 3, 0]
    @test is.three_nis == [
        Peridynamics.ThreeNeighborInteraction(1, 2, 3, 1.0),
        Peridynamics.ThreeNeighborInteraction(4, 5, 6, 1.0),
        Peridynamics.ThreeNeighborInteraction(7, 8, 9, 1.0),
        Peridynamics.ThreeNeighborInteraction(10, 11, 12, 1.0)]
    @test is.n_three_nis == [1, 1, 1, 1, 0]

    @test is.volume ≈ [1, 1, 1, 1, 1]
    @test is.volume_one_nis ≈ [fill(1.33333333333333, 4); 0.0]
    @test is.volume_two_nis ≈ [fill(1.33333333333333, 4); 0.0]
    @test is.volume_three_nis ≈ [fill(4.0, 4); 0.0]
end

@testitem "interaction system compatibility" begin
    # setup
    pos, vol = uniform_box(1, 0.25, 0.25, 0.25)
    body = Body(BBMaterial(), pos, vol)
    @test_throws ArgumentError Peridynamics.check_interaction_system_compat(body.mat)
end

@testitem "interaction counts of a uniform cube" setup=[Fixtures] begin
    # The one-, two- and three-neighbor interactions of a 4³ cube with a horizon of 2.015 Δx.
    # The counts are regression values: every point has the same family as in a bond system
    # (1128 one-neighbor interactions in total), and the pairs and triplets are built from those
    # families by `find_two_nis` and `find_three_nis`; a change here means the family or the
    # pair/triplet search changed. The interaction parameters are explicit, because with
    # `nu = 0.25` no triplets would be built at all (see `Fixtures.cki_kwargs`).
    Δx = 0.25
    pos, vol = uniform_box(1, 1, 1, Δx)
    body = Body(CKIMaterial(), pos, vol)
    @test_logs (:warn, r"specified manually") material!(body; horizon=2.015Δx, rho=7850,
                                                          E=210e9, nu=0.25,
                                                          Fixtures.cki_kwargs()...)
    pd = Peridynamics.PointDecomposition(body, 1)
    system = Peridynamics.get_system(body, pd, 1)
    @test n_points(body) == 64
    @test length(system.one_nis) == 1128
    @test length(system.two_nis) == 5400
    @test length(system.three_nis) == 9144
    # every interaction count is consistent with the per-point bookkeeping
    @test sum(system.n_one_nis) == length(system.one_nis)
    @test sum(system.n_two_nis) == length(system.two_nis)
    @test sum(system.n_three_nis) == length(system.three_nis)
end

@testitem "InteractionSystem: a precrack breaks the one-neighbor interactions" setup=[Fixtures] begin
    body = Fixtures.cube(CKIMaterial(); n=4, m=2.015)
    point_set!(p -> p[1] < 0, body, :left)
    point_set!(p -> p[1] >= 0, body, :right)
    precrack!(body, :left, :right)
    c = Fixtures.chunk(body)
    (; storage, system) = c
    @test all(storage.one_ni_active)
    Peridynamics.apply_precracks!(c, body)
    left, right = body.point_sets[:left], body.point_sets[:right]
    crossing = [(i in left) != (system.one_nis[bid].neighbor in left)
                for i in Peridynamics.each_point_idx(system)
                for bid in Peridynamics.each_one_ni_idx(system, i)]
    @test any(crossing)
    @test storage.one_ni_active == .!crossing
    # the active counts are consistent with the flags, and the damage sees the precrack
    n_active = [count(storage.one_ni_active[Peridynamics.each_one_ni_idx(system, i)])
                for i in Peridynamics.each_point_idx(system)]
    @test storage.n_active_one_nis == n_active
    @test Peridynamics.one_ni_failure(storage, findfirst(crossing)) == false
end

@testitem "InteractionSystem: required parameters and log message" begin
    params = Peridynamics.required_point_parameters(CKIMaterial)
    @test :δ in params && :rho in params && :C1 in params && :C2 in params && :C3 in params
    msg = Peridynamics.log_msg_interaction_system(11, 22, 33)
    @test contains(msg, "one-neighbor-interactions") && contains(msg, "11")
    @test contains(msg, "two-neighbor-interactions") && contains(msg, "22")
    @test contains(msg, "three-neighbor-interactions") && contains(msg, "33")
end
