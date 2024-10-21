
function mode_i(mat, npyz)
    l, Δx, δ, a = 1.0, 1/npyz, 3.015/npyz, 0.5
    pos, vol = uniform_box(l, l, 0.1l, Δx)
    ids = sortperm(pos[2,:])
    body = Body(mat, pos[:, ids], vol[ids])
    material!(body; horizon=3.015Δx, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
    point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, body, :set_a)
    point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, body, :set_b)
    precrack!(body, :set_a, :set_b)
    point_set!(p -> p[2] > l/2-Δx, body, :set_top)
    point_set!(p -> p[2] < -l/2+Δx, body, :set_bottom)
    velocity_bc!(t -> -30, body, :set_bottom, :y)
    velocity_bc!(t -> 30, body, :set_top, :y)
    vv = VelocityVerlet(steps=2000)
    job = Job(body, vv; path=joinpath(@__DIR__, "results", "mode_i-BB-npyz$(npyz)"))
    return @benchmarkable submit($job)
end
