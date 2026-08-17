"""
    wave_in_bar(; kwargs...)

A single sinusoidal velocity pulse travelling along a slender bar.

The front of the pulse moves at the one-dimensional wave speed `c₀ = sqrt(E / rho)`, so the
measured speed can be compared against a closed form without a reference solution. That makes
this the cheapest check that the physics is still right, and refining `npyz` turns it into a
convergence study. The default `time` stops before the pulse reaches the far end, so that no
reflection contaminates the measurement.

# Keywords
- `mat`: Material of the body. (default: `BBMaterial()`)
- `npyz::Int`: Number of points over the cross section. (default: `5`)
- `m::Real`: Ratio of horizon to point spacing. (default: `3.015`)
- `lx`, `lyz`: Length and edge length of the square cross section. (default: `0.2`, `0.002`)
- `E`, `nu`, `rho`: Material parameters. Note that the bond-based materials are restricted to
    `nu = 1/4`. (defaults: steel)
- `vmax::Real`: Peak velocity of the pulse. (default: `2.0`)
- `T::Real`: Duration of the pulse. (default: `1e-5`)
- `time::Real`: Simulated time. (default: `4e-5`)
- `steps`: Number of time steps. Overrides `time` if given, so that the amount of work is
    fixed instead of the simulated duration. (default: `nothing`)
- `sort_points::Bool`: Sort the points along the bar axis. Only matters with MPI, where it is
    essential: `uniform_box` varies `x` fastest, so a chunk of consecutive indices is a handful
    of complete lines along the bar and its halo is almost the whole body. Sorting by `x` turns
    the same decomposition into slabs whose halo is two thin slices. Off by default so that
    timings stay comparable with baselines measured before it existed. (default: `false`)
- `freq::Int`: Export frequency, only relevant together with `path`. (default: `10`)
- `path`: Where to export results. Nothing is exported if this is `nothing`. (default:
    `nothing`)
"""
function wave_in_bar(; mat=BBMaterial(), npyz::Int=5, m::Real=3.015, lx::Real=0.2,
                     lyz::Real=0.002, E::Real=210e9, nu::Real=0.25, rho::Real=7850.0,
                     vmax::Real=2.0, T::Real=1e-5, time::Real=4e-5, steps=nothing,
                     sort_points::Bool=false, freq::Int=10, path=nothing)
    Δx = lyz / npyz
    pos, vol = uniform_box(lx, lyz, lyz, Δx)
    if sort_points
        ids = sortperm(@view pos[1, :])
        pos, vol = pos[:, ids], vol[ids]
    end
    body = Body(mat, pos, vol)
    material!(body; horizon=m * Δx, rho, E, nu)
    no_failure!(body)
    point_set!(p -> p[1] < -lx / 2 + 1.2Δx, body, :left)
    velocity_bc!(t -> t < T ? vmax * sinpi(2t / T) : 0.0, body, :left, :x)
    solver = isnothing(steps) ? VelocityVerlet(; time, safety_factor=0.7) :
             VelocityVerlet(; steps, safety_factor=0.7)
    isnothing(path) && return Job(body, solver)
    return Job(body, solver; path, freq, fields=(:displacement, :velocity))
end
