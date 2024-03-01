let
    l, Δx, δ, a = 1.0, 1/4, 3.015/4, 0.5
    pos, vol = uniform_box(l, l, l, Δx)
    b = Body(BBMaterial(), pos, vol)
    material!(b; horizon=3.015Δx, E=2.1e5, rho=8e-6, Gc=2.7)
    point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
    point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
    precrack!(b, :set_a, :set_b)
    point_set!(p -> p[2] > l/2-Δx, b, :set_top)
    point_set!(p -> p[2] < -l/2+Δx, b, :set_bottom)
    velocity_ic!(b, :set_bottom, :x, 0)
    velocity_ic!(b, :set_bottom, :y, 0)
    velocity_ic!(b, :set_bottom, :z, 0)
    velocity_ic!(b, :set_top, 1, 0)
    velocity_ic!(b, :set_top, 2, 0)
    velocity_ic!(b, :set_top, 3, 0)
    velocity_bc!(t -> 0, b, :set_bottom, :x)
    velocity_bc!(t -> -30, b, :set_bottom, :y)
    velocity_bc!(t -> 0, b, :set_bottom, :z)
    velocity_bc!(t -> 0, b, :set_top, 1)
    velocity_bc!(t -> 30, b, :set_top, 2)
    velocity_bc!(t -> 0, b, :set_top, 3)
    forcedensity_bc!(t -> 0, b, :set_top, :x)
    forcedensity_bc!(t -> 0, b, :set_top, :y)
    forcedensity_bc!(t -> 0, b, :set_top, :z)
    forcedensity_bc!(t -> 0, b, :set_bottom, 1)
    forcedensity_bc!(t -> 0, b, :set_bottom, 2)
    forcedensity_bc!(t -> 0, b, :set_bottom, 3)
    vv = VelocityVerlet(steps=1)
    job = Job(b, vv)
    submit(job; quiet=true);
end

let
    l, Δx, δ, a = 1.0, 1/4, 3.015/4, 0.5
    pos, vol = uniform_box(l, l, l, Δx)
    b = Body(OSBMaterial(), pos, vol)
    material!(b; horizon=3.015Δx, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
    point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
    point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
    failure_permit!(b, :set_a, false)
    precrack!(b, :set_a, :set_b)
    point_set!(p -> p[2] > l/2-Δx, b, :set_top)
    point_set!(p -> p[2] < -l/2+Δx, b, :set_bottom)
    velocity_bc!(t -> -30, b, :set_bottom, :y)
    velocity_bc!(t -> 30, b, :set_top, 2)
    vv = VelocityVerlet(steps=1)
    job = Job(b, vv)
    submit(job; quiet=true);
end

let
    l, Δx, δ, a = 1.0, 1/4, 3.015/4, 0.5
    pos, vol = uniform_box(l, l, l, Δx)
    b = Body(NOSBMaterial(), pos, vol)
    material!(b; horizon=3.015Δx, E=2.1e5, nu=0.25, rho=8e-6, Gc=2.7)
    point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
    point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
    precrack!(b, :set_a, :set_b)
    point_set!(p -> p[2] > l/2-Δx, b, :set_top)
    point_set!(p -> p[2] < -l/2+Δx, b, :set_bottom)
    velocity_bc!(t -> -30, b, :set_bottom, :y)
    velocity_bc!(t -> 30, b, :set_top, 2)
    vv = VelocityVerlet(steps=1)
    job = Job(b, vv)
    submit(job; quiet=true);
end
