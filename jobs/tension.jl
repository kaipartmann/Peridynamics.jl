"""
    tension(; kwargs...)

A uniform cube pulled apart along the x-axis.

The kernel benchmarks and the zero-energy mode verification are built on this body, because a
uniform cube gives every interior point a full neighborhood, so the answer away from the
surfaces is the classical one and can be compared against a closed form. Both impose the
deformation directly on the chunk instead of using the condition below, because that prescribes
the field rather than solving for it: conditions are applied by the time solver and not by the
force density kernel, which is the only thing those two look at.

# Keywords
- `mat`: Material of the body. (default: `BBMaterial()`)
- `npyz::Int`: Number of points along each edge. (default: `8`)
- `m::Real`: Ratio of horizon to point spacing. (default: `3.015`)
- `E`, `nu`, `rho`, `Gc`: Material parameters. Note that the bond-based materials are
    restricted to `nu = 1/4`. (defaults: steel)
- `v::Real`: Velocity imposed on each of the two end faces, in opposite directions.
    (default: `1.0`)
- `failure::Bool`: Whether bonds may break. Off by default, because the critical stretch is of
    the order of `1e-5` here, so any deformation worth measuring breaks every bond at once.
    Note that `calc_failure!` still visits every bond either way, so the cost stays
    representative. (default: `false`)
- `steps::Int`: Number of time steps. (default: `1`)
- `path`: Where to export results. Nothing is exported if this is `nothing`. (default:
    `nothing`)
"""
function tension(; mat=BBMaterial(), npyz::Int=8, m::Real=3.015, E::Real=210e9,
                 nu::Real=0.25, rho::Real=7850.0, Gc::Real=2.7, v::Real=1.0,
                 failure::Bool=false, steps::Int=1, path=nothing)
    l = 1.0
    Δx = l / npyz
    pos, vol = uniform_box(l, l, l, Δx)
    body = Body(mat, pos, vol)
    material!(body; horizon=m * Δx, rho, E, nu, Gc)
    failure || no_failure!(body)
    point_set!(p -> p[1] > l / 2 - Δx, body, :right)
    point_set!(p -> p[1] < -l / 2 + Δx, body, :left)
    velocity_bc!(t -> v, body, :right, :x)
    velocity_bc!(t -> -v, body, :left, :x)
    solver = VelocityVerlet(; steps)
    isnothing(path) && return Job(body, solver)
    return Job(body, solver; path)
end
