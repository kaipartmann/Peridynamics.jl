@testitem "one- and two-neighbor interactions" begin
    position = [0.0 1.0 0.0
                0.0 0.0 0.0
                0.0 0.0 1.0]
    volume = fill(1.0, 3)
    body = Body(CKIMaterial(), position, volume)
    material!(body; horizon=1.5, rho=8e-6, E=210e3, nu=0.3, Gc=1.0, C1=1, C2=1)
    pd = Peridynamics.PointDecomposition(body, 1)

    is, ch = Peridynamics.InteractionSystem(body, pd, 1)

    @test is.one_nis == [
        Peridynamics.Bond(2, 1.0, true),
        Peridynamics.Bond(3, 1.0, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(3, √2, true),
        Peridynamics.Bond(1, 1.0, true),
        Peridynamics.Bond(2, √2, true)]
    @test is.n_one_nis == [2, 2, 2]
    @test is.one_ni_idxs == [1:2, 3:4, 5:6]
    @test is.two_nis == [
        Peridynamics.TwoNeighborInteraction(1, 2, 1.0),
        Peridynamics.TwoNeighborInteraction(3, 4, 1.0),
        Peridynamics.TwoNeighborInteraction(5, 6, 1.0)]
    @test is.n_two_nis == [1, 1, 1]
    @test is.two_ni_idxs == [1:1, 2:2, 3:3]
    @test is.three_nis == Vector{Peridynamics.ThreeNeighborInteraction}()
    @test is.n_three_nis == Vector{Int}()
    @test is.three_ni_idxs == Vector{UnitRange{Int}}()
end

@testitem "one-, two-, and three-neighbor interactions" begin
    position = [0.0 1.0 0.0 0.0 2.0
                0.0 0.0 1.0 0.0 2.0
                0.0 0.0 0.0 1.0 2.0]
    volume = fill(1.0, 5)
    body = Body(CKIMaterial(), position, volume)
    material!(body; horizon=1.5, rho=8e-6, E=210e3, nu=0.3, Gc=1.0, C1=1, C2=1, C3=1)
    pd = Peridynamics.PointDecomposition(body, 1)

    is, ch = Peridynamics.InteractionSystem(body, pd, 1)

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
    @test is.two_ni_idxs == [1:3, 4:6, 7:9, 10:12, 0:-1]
    @test is.three_nis == [
        Peridynamics.ThreeNeighborInteraction(1, 2, 3, 1.0),
        Peridynamics.ThreeNeighborInteraction(4, 5, 6, 1.0),
        Peridynamics.ThreeNeighborInteraction(7, 8, 9, 1.0),
        Peridynamics.ThreeNeighborInteraction(10, 11, 12, 1.0)]
    @test is.n_three_nis == [1, 1, 1, 1, 0]
    @test is.three_ni_idxs == [1:1, 2:2, 3:3, 4:4, 0:-1]
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
