"""
    torsional_frequency(E, nu, rho, l)

Angular frequency of the first torsional mode of a clamped-free rod of length `l`, from the
classical wave equation `θ_tt = c² θ_xx` with `c = sqrt(G / rho)`.

The cross section does not enter: for the fundamental mode it cancels between the polar moment
of inertia and the torsional constant, which is what makes this the reference for a peridynamic
bar of any thickness.
"""
function torsional_frequency(E::Real, nu::Real, rho::Real, l::Real)
    G = E / (2 * (1 + nu))
    return (π / 2) * sqrt(G / rho) / l
end

"""
    bartwist(; kwargs...)

A slender bar, clamped at one end, released from the velocity field of its first torsional mode
and left to swing freely about its own axis.

This is the large rotation verification case. The free end turns by more than a radian, far
outside the small rotation regime, while the strains along the way stay of the order of a
percent: a formulation that is not objective mistakes part of that rotation for deformation,
which shows up as a torsional frequency that no longer matches [`torsional_frequency`](@ref).
Correspondence formulations are the ones at risk here, they build their stress from a
deformation gradient rather than from bond stretches.

The bar is aligned with the x-axis. The clamped part is added on top of `lz` and is at least one
horizon thick, so that the free length is `lz` and the clamp is not a single layer of points
that the rest of the bar reaches across.

# Keywords
- `mat`: Material of the body. (default: `BBMaterial()`)
- `npxy::Int`: Number of points across the square cross section. (default: `10`)
- `lxy::Real`: Edge length of the square cross section. (default: `0.1`)
- `lz::Real`: Free length of the bar, rounded to a whole number of point spacings. Whoever
    compares against [`torsional_frequency`](@ref) has to use the same length, so pick a `lz`
    that is a multiple of `lxy / npxy`. (default: `0.5`)
- `m::Real`: Ratio of horizon to point spacing. (default: `3.015`)
- `E`, `nu`, `rho`, `Gc`: Material parameters. Note that the bond-based materials are
    restricted to `nu = 1/4`. The default is a soft polymer, which keeps the torsional wave
    speed and therefore the number of time steps low. `Gc` only exists because the materials
    with a damage model ask for it, bonds never break here.
- `Ω_0::Real`: Peak angular velocity of the initial mode shape, reached at the free end. The
    default turns the free end by about 1.2 rad, or 70 degrees. (default: `300.0`)
- `time::Real`: Simulated time. The default is one analytical period, which is enough to see
    the free end swing out and come back through zero even if the discretization softens the
    bar considerably. (default: one period)
- `freq::Int`: Export frequency, only relevant together with `path`. (default: `10`)
- `path`: Where to export results. Nothing is exported if this is `nothing`. (default:
    `nothing`)
"""
function bartwist(; mat=BBMaterial(), npxy::Int=10, lxy::Real=0.1, lz::Real=0.5,
                  m::Real=3.015, E::Real=17e6, nu::Real=0.25, rho::Real=1100.0,
                  Gc::Real=100.0, Ω_0::Real=300.0,
                  time::Real=2π / torsional_frequency(E, nu, rho, lz), freq::Int=10,
                  path=nothing)
    Δx = lxy / npxy
    n_clamp = ceil(Int, m)
    n_free = round(Int, lz / Δx)
    # the free length is what the analytical frequency is taken over, so it has to be the one
    # the grid actually realizes and not the requested `lz`
    l = n_free * Δx
    lx = (n_free + n_clamp) * Δx
    pos, vol = uniform_box(lx, lxy, lxy, Δx; center=(-n_clamp * Δx / 2, 0, 0))
    body = Body(mat, pos, vol)
    material!(body; horizon=m * Δx, E, nu, rho, Gc)
    # the bar is twisted well beyond any critical stretch of a real material, and a broken bond
    # would end the vibration this case is about
    no_failure!(body)
    point_set!(p -> p[1] < -l / 2, body, :clamped)
    point_set!(p -> p[1] ≥ -l / 2, body, :free)
    for dim in (:x, :y, :z)
        velocity_bc!(t -> 0.0, body, :clamped, dim)
    end
    # rigid rotation about the x-axis at the rate `Ω_0 * sin(π x / 2l)`, the first mode shape
    velocity_ic!(body, :free, :y) do p
        return -Ω_0 * sinpi((p[1] + l / 2) / 2l) * p[3]
    end
    velocity_ic!(body, :free, :z) do p
        return +Ω_0 * sinpi((p[1] + l / 2) / 2l) * p[2]
    end
    solver = VelocityVerlet(; time)
    isnothing(path) && return Job(body, solver)
    return Job(body, solver; path, freq, fields=(:displacement,))
end
