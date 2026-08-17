"""
    mode_i(; kwargs...)

A pre-cracked plate pulled apart in mode I, in which a single straight crack runs from the notch
to the far edge.

This is the whole-simulation benchmark and the damage model verification case: it exercises the
failure bookkeeping and the time integration together with the force density kernel, and the
speed of its crack is a material property that can be compared against the Rayleigh wave speed.

The parameters keep the applied far-field strain below the critical stretch for the whole
simulation, so that only the notch tip reaches it. That is what keeps the crack clean; raising
`v` or `Gc` far from the defaults instead breaks bonds along the driven edges, and the plate
then fails at its boundary.

# Keywords
- `mat`: Material of the body. (default: `BBMaterial()`)
- `npxy::Int`: Number of points along the plate edge. The plate is `0.1l` thick.
    (default: `40`)
- `l::Real`: Edge length of the square plate. (default: `0.1`)
- `a::Real`: Length of the precrack, measured from the left edge. (default: `0.05`)
- `m::Real`: Ratio of horizon to point spacing. (default: `3.015`)
- `E`, `nu`, `rho`, `Gc`: Material parameters. Note that the bond-based materials are
    restricted to `nu = 1/4`.
- `v::Real`: Velocity imposed on the top and bottom edges. (default: `0.1`)
- `time::Real`: Simulated time. The default is long enough for the crack to cross the plate
    and arrest at the far edge. (default: `3e-4`)
- `freq::Int`: Export frequency, only relevant together with `path`. (default: `10`)
- `path`: Where to export results. Nothing is exported if this is `nothing`. (default:
    `nothing`)
"""
function mode_i(; mat=BBMaterial(), npxy::Int=40, l::Real=0.1, a::Real=0.05, m::Real=3.015,
                E::Real=32e9, nu::Real=0.25, rho::Real=2500.0, Gc::Real=100.0, v::Real=0.1,
                time::Real=3e-4, freq::Int=10, path=nothing)
    Δx = l / npxy
    δ = m * Δx
    pos, vol = uniform_box(l, l, 0.1l, Δx)
    ids = sortperm(pos[2, :])
    body = Body(mat, pos[:, ids], vol[ids])
    material!(body; horizon=δ, E, nu, rho, Gc)
    point_set!(p -> p[1] ≤ -l / 2 + a && 0 ≤ p[2] ≤ 2δ, body, :set_a)
    point_set!(p -> p[1] ≤ -l / 2 + a && -2δ ≤ p[2] < 0, body, :set_b)
    precrack!(body, :set_a, :set_b)
    point_set!(p -> p[2] > l / 2 - Δx, body, :set_top)
    point_set!(p -> p[2] < -l / 2 + Δx, body, :set_bottom)
    velocity_bc!(t -> -v, body, :set_bottom, :y)
    velocity_bc!(t -> v, body, :set_top, :y)
    solver = VelocityVerlet(; time)
    isnothing(path) && return Job(body, solver)
    return Job(body, solver; path, freq, fields=(:displacement, :damage))
end
