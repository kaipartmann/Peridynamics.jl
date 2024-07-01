@testitem "Point duplicates" begin
    pos, vol = uniform_box(1, 1, 1, 1/5)
    pos[:, end] .= pos[:, end-1]
    body = Body(BBMaterial(), pos, vol)
    material!(body, horizon=4/100, E=1, rho=1, Gc=1)
    @test_throws ErrorException Peridynamics.find_bonds(body, 1:n_points(body))
end
