# The simulations of the MPI-vs-threads comparison. This plain file is included both by the
# `MPIComparison` test module (for the threads run in the test process) and by `comparison.jl`
# (for the MPI run in the `mpiexec` subprocess), so both sides cannot drift apart.
#
# Each case is a small plate with a pre-crack, pulled apart at the top and the bottom. The
# decomposition of a simulation into chunks must not change the result, so the threads run and
# the MPI run are compared file by file. That is a property of the decomposition and the halo
# exchange, not of the time integration, so the runs are kept short: the dynamic ones stop after
# 20 steps, which is enough to involve every chunk and every halo exchange path.

"Points per edge of the plate; the plate has `N × N × N/10` points."
const COMPARISON_N = 15

"""
    comparison_body(mat, l; precrack=true, kwargs...)

A plate of size `l × l × l/10`, sorted along `y` so that the decomposition splits it into
horizontal stripes, with a pre-crack from the left edge to the center unless `precrack=false`.
`kwargs` are passed to `material!`.
"""
function comparison_body(mat, l; precrack=true, kwargs...)
    Δx, δ, a = l / COMPARISON_N, 3.015 * l / COMPARISON_N, 0.5l
    pos, vol = uniform_box(l, l, 0.1l, Δx)
    ids = sortperm(pos[2, :])
    body = Body(mat, pos[:, ids], vol[ids])
    material!(body; horizon=3.015Δx, kwargs...)
    if precrack
        point_set!(p -> p[1] ≤ -l/2 + a && 0 ≤ p[2] ≤ 2δ, body, :set_a)
        point_set!(p -> p[1] ≤ -l/2 + a && -2δ ≤ p[2] < 0, body, :set_b)
        precrack!(body, :set_a, :set_b)
    end
    point_set!(p -> p[2] > l/2 - Δx, body, :set_top)
    point_set!(p -> p[2] < -l/2 + Δx, body, :set_bottom)
    return body
end

"A `VelocityVerlet` pull-apart of `mat`, 20 steps, exported every 10 steps (3 files)."
function dynamic_case(mat; kwargs...)
    body = comparison_body(mat, 1.0; E=2.1e5, rho=8e-6, Gc=2.7, kwargs...)
    velocity_bc!(t -> -30, body, :set_bottom, :y)
    velocity_bc!(t -> 30, body, :set_top, :y)
    return body, VelocityVerlet(steps=20), 10
end

"""
A `DynamicRelaxation` pull-apart, 200 steps, exported at the end (2 files).

No pre-crack here: `calc_damping` estimates the damping coefficient per chunk from the local
degrees of freedom, so the relaxation is only independent of the decomposition when every chunk
sees the same deformation state, i.e. for a uniform field. A crack makes the chunks differ and
the threads and MPI results with them; that is a property of the per-chunk damping estimate,
checked in the verification suite ("decomposition invariance"), and not of the halo exchange.
"""
function relaxation_case(mat; kwargs...)
    body = comparison_body(mat, 100.0; precrack=false, E=2.1e5, rho=8e-6, kwargs...)
    forcedensity_bc!(t -> -1e3, body, :set_bottom, :y)
    forcedensity_bc!(t -> 1e3, body, :set_top, :y)
    return body, DynamicRelaxation(steps=200, damping_factor=0.05), 200
end

"A `NewtonKrylov` pull-apart, 5 load steps, exported at the end (2 files)."
function newton_case(mat; kwargs...)
    body = comparison_body(mat, 100.0; E=2.1e5, rho=8e-6, kwargs...)
    forcedensity_bc!(p -> -1e3, body, :set_bottom, :y)
    forcedensity_bc!(p -> 1e3, body, :set_top, :y)
    return body, NewtonKrylov(steps=5, stepsize=1.0, maxiter=200, tol=1e-3), 5
end

"""
Name => (case builder, number of exported files, comparison tolerance). Every material that has
an MPI-relevant code path of its own (corrections, gradient weights, bond-associated families)
is here. `nothing` as tolerance means the results must agree to `≈`.
"""
const COMPARISON_CASES = [
    "BB VelocityVerlet" => (() -> dynamic_case(BBMaterial()), 3, nothing),
    "BB DynamicRelaxation" => (() -> relaxation_case(BBMaterial()), 2, nothing),
    # The Newton-Krylov comparison is loose on purpose, see
    # https://github.com/kaipartmann/Peridynamics.jl/issues/187
    "BB NewtonKrylov" => (() -> newton_case(BBMaterial()), 2, 0.03),
    "BB-ESC VelocityVerlet" => (() -> dynamic_case(BBMaterial{EnergySurfaceCorrection}()), 3,
                                nothing),
    "OSB VelocityVerlet" => (() -> dynamic_case(OSBMaterial(); nu=0.25), 3, nothing),
    "C VelocityVerlet" => (() -> dynamic_case(CMaterial(); nu=0.25), 3, nothing),
    "RKCR VelocityVerlet" => (() -> dynamic_case(RKCRMaterial(); nu=0.25), 3, nothing),
    "BAC VelocityVerlet" => (() -> dynamic_case(BACMaterial(); nu=0.25), 3, nothing),
]

"""
    run_comparison_case(name, path; n_chunks=nothing)

Run the case `name` and export it to `path`. Under `mpiexec` every rank owns one chunk; with
threads the body is split into `n_chunks` chunks (the number of threads if `nothing`). The
dynamic relaxation is only independent of the decomposition if the chunks see the same
deformation state (see `relaxation_case`), so the comparison runs both sides with
`COMPARISON_NRANKS` chunks: a symmetric top/bottom split of the plate.
"""
function run_comparison_case(name::AbstractString, path::AbstractString; n_chunks=nothing)
    builder, _, _ = Dict(COMPARISON_CASES)[name]
    body, solver, freq = builder()
    job = Job(body, solver; path, freq)
    if n_chunks === nothing
        submit(job; quiet=true)
    else
        Peridynamics.set_quiet!(true)
        Peridynamics.submit_threads(job, n_chunks)
    end
    return nothing
end

"Number of MPI ranks (and threads chunks) of the comparison."
const COMPARISON_NRANKS = 4

"The directory the results of case `name` are exported to below `root`."
case_path(root, name) = joinpath(root, replace(name, " " => "_"))
