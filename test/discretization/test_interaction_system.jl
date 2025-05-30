@testitem "one- and two-neighbor interactions" begin
    position = [0.0 1.0 0.0
                0.0 0.0 0.0
                0.0 0.0 1.0]
    volume = fill(1.0, 3)
    body = Body(CKIMaterial(), position, volume)
    material!(body; horizon=1.5, rho=8e-6, E=210e3, nu=0.3, Gc=1.0, C1=1, C2=1)
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
    material!(body, :a; horizon=1.5, rho=7e-6, E=200e3, nu=0.3, Gc=1.0, C1=1, C2=1)

    ts = VelocityVerlet(steps=1) # nur als dummy
    dh = Peridynamics.threads_data_handler(body, ts, 1) # ein thread
    chunk = dh.chunks[1] # das erste und einzige chunk
    @test Peridynamics.has_three_nis(chunk.paramsetup) == false
end

@testitem "one- and three-neighbor interactions" begin
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 0.0 1.0
                0.0 0.0 1.0 0.0]
    volume = fill(1.0, 4)
    body = Body(CKIMaterial(), position, volume)
    material!(body; horizon=1.5, rho=8e-6, E=210e3, nu=0.3, Gc=1.0, C1=1, C3=1)
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
    material!(body, :a; horizon=1.5, rho=7e-6, E=200e3, nu=0.3, Gc=1.0, C1=1, C3=1)

    ts = VelocityVerlet(steps=1) # nur als dummy
    dh = Peridynamics.threads_data_handler(body, ts, 1) # ein thread
    chunk = dh.chunks[1] # das erste und einzige chunk
    @test Peridynamics.has_two_nis(chunk.paramsetup) == false
end

@testitem "one-, two-, and three-neighbor interactions" begin
    position = [0.0 1.0 0.0 0.0 2.0
                0.0 0.0 1.0 0.0 2.0
                0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    body = Body(CKIMaterial(), position, volume)
    material!(body; horizon=1.5, rho=8e-6, E=210e3, nu=0.3, Gc=1.0, C1=1, C2=1, C3=1)
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
    material!(body; horizon=1.5, rho=8e-6, E=210e3, nu=0.3, Gc=1.0, C1=1, C2=1, C3=1)

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
