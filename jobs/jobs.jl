# Job definitions shared by `benchmark/` and `test/verification/`, so that both of them measure
# the same simulation. Every file here defines one function that returns a `Job`, with all
# arguments as keywords. Load it with
#
#     include(joinpath(pkgdir(Peridynamics), "jobs", "jobs.jl"))
#
# Nothing here may depend on anything but Peridynamics and the standard library, and nothing
# here may use the storage or parameter macros. Both keep this directory working across changes
# to the internals.

using Peridynamics

include(joinpath(@__DIR__, "tension.jl"))
include(joinpath(@__DIR__, "wave_in_bar.jl"))
include(joinpath(@__DIR__, "mode_i.jl"))

"""
    ChunkFixture

A body chunk with the time and step size at which its force density is evaluated. The step size
is part of it because the rotated formulations update their stress history incrementally, so an
arbitrary value changes both the state and the code path.
"""
struct ChunkFixture{C}
    chunk::C
    t::Float64
    Δt::Float64
end

"""
    chunk_fixture(job)

The single initialized body chunk of `job` and its stable time step, exactly as a time solver
would see them, but with no threading, MPI, file output or logging around it.
"""
function chunk_fixture(job)
    dh = Peridynamics.threads_data_handler(job.spatial_setup, job.time_solver, 1)
    Peridynamics.init_time_solver!(job.time_solver, dh)
    Peridynamics.initialize!(dh, job.time_solver)
    Δt = job.time_solver.Δt
    return ChunkFixture(dh.chunks[1], Δt, Δt)
end

"""
    deform!(chunk; F, amplitude, Δt)

Impose a small non-identity `F` plus an inhomogeneous perturbation on `chunk` and return it.

Both parts matter: an undeformed body has zero bond strain and short-circuits, and a purely
homogeneous deformation makes the zero-energy mode stabilization of the correspondence
formulations vanish, which is the expensive part of exactly those materials. The perturbation is
deterministic so that timings stay comparable between runs.

With `Δt`, the displacement is treated as an increment over that step and the velocity fields are
populated too, which `CRMaterial` and `RKCRMaterial` need for the velocity gradient.
"""
function deform!(chunk; F=[1.01 0.003 0.0; 0.0 0.995 0.0; 0.0 0.0 1.0], amplitude=1e-4,
                 Δt=nothing)
    ref = chunk.system.position
    (; position, displacement, velocity, velocity_half) = chunk.storage
    for i in axes(ref, 2)
        X1, X2, X3 = ref[1, i], ref[2, i], ref[3, i]
        x1 = F[1, 1] * X1 + F[1, 2] * X2 + F[1, 3] * X3 + amplitude * sinpi(3X2)
        x2 = F[2, 1] * X1 + F[2, 2] * X2 + F[2, 3] * X3 + amplitude * sinpi(3X3)
        x3 = F[3, 1] * X1 + F[3, 2] * X2 + F[3, 3] * X3 + amplitude * sinpi(3X1)
        for (d, x, X) in ((1, x1, X1), (2, x2, X2), (3, x3, X3))
            u = x - X
            position[d, i] = x
            displacement[d, i] = u
            if !isnothing(Δt)
                velocity[d, i] = u / Δt
                velocity_half[d, i] = u / Δt
            end
        end
    end
    return chunk
end

function deform!(fixture::ChunkFixture; kwargs...)
    deform!(fixture.chunk; kwargs..., Δt=fixture.Δt)
    return fixture
end

"""
    force_density!(fixture)

Run the complete force density calculation of `fixture`, the hot loop of every simulation.

The reproducing kernel materials override `calc_force_density!` for the data handler and not
for the chunk, because they need a halo exchange between the gradient weights and the forces.
Calling the chunk method alone would silently skip `calc_weights_and_defgrad!` and measure half
the work, so those materials get their own method here.
"""
force_density!(fixture::ChunkFixture) = force_density!(fixture.chunk, fixture.t, fixture.Δt)

force_density!(chunk, t, Δt) = Peridynamics.calc_force_density!(chunk, t, Δt)

function force_density!(chunk::Peridynamics.BodyChunk{<:Peridynamics.AbstractBondSystem,
                                                      <:Peridynamics.AbstractRKCMaterial},
                        t, Δt)
    Peridynamics.calc_weights_and_defgrad!(chunk, t, Δt)
    Peridynamics.calc_force_density!(chunk, t, Δt)
    return nothing
end

"""
    gradient_weights!(chunk)

Recompute the gradient weights of a reproducing kernel material for every point of `chunk`.

[`force_density!`](@ref) does not reach this. The weights are only recomputed where damage has
just grown, and `calc_damage!` rewrites the `update_gradients` flag at the start of every force
calculation, so an undamaged body never enters it and the cost stays invisible.
"""
gradient_weights!(chunk) = Peridynamics.initialize!(chunk)
gradient_weights!(fixture::ChunkFixture) = gradient_weights!(fixture.chunk)

"Whether [`gradient_weights!`](@ref) measures anything for `mat`."
has_gradient_weights(mat) = mat isa Peridynamics.AbstractRKCMaterial

# Resolution of the body a material is measured on. The interaction system materials evaluate
# triplets of neighbors instead of pairs, so they need a smaller body and a smaller horizon to
# stay in the same order of magnitude as the rest.
fixture_size(::Peridynamics.AbstractMaterial) = (npyz=8, m=3.015)
fixture_size(::Peridynamics.AbstractInteractionSystemMaterial) = (npyz=5, m=2.015)

"A deformed [`tension`](@ref) chunk fixture for `mat`, ready to be measured."
material_fixture(mat) = deform!(chunk_fixture(tension(; mat, fixture_size(mat)...)))

"""
    rayleigh_wave_speed(E, nu, rho)

Rayleigh wave speed of an isotropic solid, from the approximation of Viktorov (1967). A crack
in that solid cannot exceed it, so it is the reference every measured crack speed is expressed
against.
"""
function rayleigh_wave_speed(E::Real, nu::Real, rho::Real)
    c_s = sqrt(E / (2 * (1 + nu)) / rho)
    return c_s * (0.862 + 1.14nu) / (1 + nu)
end

"""
    crack_tip(position, damage; dmg_min, dmg_max, xmin, cluster_radius)

Position of the leading crack tip as `(x, |y|)`, or `(NaN, NaN)` if there is no crack.

`damage` is taken in a band and not above a threshold, because a crack tip is *partially*
damaged: behind it the material is separated and ahead of it there is diffuse micro damage that
is not a crack. `xmin` keeps the search ahead of the precrack, whose points are fully damaged
from the first step, and averaging over a cluster keeps a single stray point from defining the
tip. The second return value is the largest `|y|` in that cluster, which says whether the crack
is still running on the symmetry line.
"""
function crack_tip(position, damage; dmg_min::Real=0.35, dmg_max::Real=0.9, xmin::Real=0.0,
                   cluster_radius::Real)
    ids = [i for i in axes(position, 2)
           if position[1, i] ≥ xmin && dmg_min ≤ damage[i] ≤ dmg_max]
    isempty(ids) && return (NaN, NaN)
    lead = ids[argmax(position[1, i] for i in ids)]
    cluster = [i for i in ids
               if sqrt(sum((position[d, i] - position[d, lead])^2 for d in 1:3)) <
                  cluster_radius]
    x = sum(position[1, i] for i in cluster) / length(cluster)
    y = maximum(abs(position[2, i]) for i in cluster)
    return (x, y)
end

"""
    crack_tip_history(job; kwargs...)

Track the crack tip of `job` over time, returning `(times, tip_x, tip_y)` sorted by time.

`job` needs a `path`, because the history is read back from its export files. That is one
simulation for the whole history, against one per sample if the job were re-run to a series of
end times.
"""
function crack_tip_history(job; cluster_radius, kwargs...)
    submit(job; quiet=true)
    times, tip_x, tip_y = Float64[], Float64[], Float64[]
    process_each_export(job, nothing; serial=true) do _, r, _
        x, y = crack_tip(r[:position], r[:damage]; cluster_radius, kwargs...)
        push!(times, first(r[:time]))
        push!(tip_x, x)
        push!(tip_y, y)
        return nothing
    end
    p = sortperm(times)
    return (times[p], tip_x[p], tip_y[p])
end

"""
    crack_speed(times, tip_x, gauge_1, gauge_2)

Speed of the crack tip between the gauge positions `gauge_1` and `gauge_2`.

Taking the speed over a fixed distance is what makes this a material property. The history has
three phases — the crack stands still while the notch loads up, then it runs, then it arrests at
the far edge — and averaging over all of them mostly measures how long the first and the last
one lasted: for the default `mode_i` that gives 0.12 c_R instead of the 0.54 c_R the crack
reaches. Both gauges therefore have to sit inside the running phase.
"""
function crack_speed(times, tip_x, gauge_1::Real, gauge_2::Real)
    t1 = gauge_crossing(times, tip_x, gauge_1)
    t2 = gauge_crossing(times, tip_x, gauge_2)
    return (gauge_2 - gauge_1) / (t2 - t1)
end

"Time at which `tip_x` first reaches `gauge`, linearly interpolated between two samples."
function gauge_crossing(times, tip_x, gauge::Real)
    i = findfirst(x -> isfinite(x) && x ≥ gauge, tip_x)
    # the samples before the crack exists are `NaN`, and the crossing cannot be interpolated
    # without a sample on the other side of the gauge
    (isnothing(i) || i == 1 || !isfinite(tip_x[i - 1])) && return NaN
    return times[i - 1] +
           (gauge - tip_x[i - 1]) / (tip_x[i] - tip_x[i - 1]) * (times[i] - times[i - 1])
end
