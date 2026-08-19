# Bodies, chunks, data handlers and deformations. Everything here returns fresh objects, because
# test items mutate what they get.

#===== bodies =====#

"""
    tetra4(mat=BBMaterial(); horizon=2.0, rho=1.0, E=1.0, nu=0.25, kwargs...)

The smallest body with a non-trivial neighborhood: the origin and the three unit vectors, with
volumes `[1.1, 1.2, 1.3, 1.4]`. Every point sees every other point with the default horizon.
`kwargs` are passed on to `material!`.
"""
function tetra4(mat=BBMaterial(); horizon=2.0, rho=1.0, E=1.0, nu=0.25, kwargs...)
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    body = Body(mat, position, volume)
    material!(body; horizon, rho, E, nu, kwargs...)
    return body
end

"Positions of `n` points on the x-axis at `0, 1, …, n-1`."
function line_position(n::Int)
    position = zeros(3, n)
    position[1, :] .= 0:(n - 1)
    return position
end

"""
    line10(mat=BBMaterial(); horizon=1.5, kwargs...)

Ten unit-volume points on the x-axis at `0, 1, …, 9`. A deterministic stand-in for the random
ten-point bodies that the `Study`, `Job` and chunk tests used to build.
"""
function line10(mat=BBMaterial(); horizon=1.5, rho=1.0, E=1.0, nu=0.25, kwargs...)
    body = Body(mat, line_position(10), ones(10))
    material!(body; horizon, rho, E, nu, kwargs...)
    return body
end

"""
    cube(mat=BBMaterial(); n=6, m=horizon_ratio(mat), kwargs...)

A unit cube of `n³` points from `uniform_box` with horizon `m * Δx`, steel-like parameters and
`no_failure!`. The critical stretch of the default parameters is of the order of 1e-5, so any
deformation worth measuring would break every bond at once; `calc_failure!` still visits every
bond, so the code path stays representative. `kwargs` override the `material!` defaults.
"""
function cube(mat=BBMaterial(); n=6, m=horizon_ratio(mat), rho=7850.0, E=210e9, nu=0.25,
              Gc=2.7, kwargs...)
    Δx = 1.0 / n
    pos, vol = uniform_box(1.0, 1.0, 1.0, Δx)
    body = Body(mat, pos, vol)
    material!(body; horizon=m * Δx, rho, E, nu, Gc, kwargs...)
    no_failure!(body)
    return body
end

#===== chunks and data handlers =====#

"""
    chunk(body, solver=VelocityVerlet(steps=1); n_chunks=1, chunk_id=1)

The cheapest `BodyChunk` there is: no data handler, no halo exchange, no threading, no
initialization of the time solver. Use it to test systems, storages and force densities of a
single chunk. `n_chunks > 1` decomposes the body and returns chunk `chunk_id`.
"""
function chunk(body, solver=VelocityVerlet(steps=1); n_chunks=1, chunk_id=1)
    pd = Peridynamics.PointDecomposition(body, n_chunks)
    ps = Peridynamics.get_param_spec(body)
    return Peridynamics.BodyChunk(body, solver, pd, chunk_id, ps)
end

"""
    handler(body, solver=VelocityVerlet(steps=1); n_chunks=1, init=true)

A threads data handler for `body`, with the time solver initialized and `initialize!` run (the
surface corrections and gradient weights are computed), so that its chunks are in the state a
time solver sees them in at the first step. `init=false` skips both initialization steps.
"""
function handler(body, solver=VelocityVerlet(steps=1); n_chunks=1, init=true)
    dh = Peridynamics.threads_data_handler(body, solver, n_chunks)
    if init
        Peridynamics.init_time_solver!(solver, dh)
        Peridynamics.initialize!(dh, solver)
    end
    return dh
end

"""
    material_fixture(mat; n=6, m=horizon_ratio(mat), kwargs...)

A deformed body chunk for `mat` together with the time and the step size at which its force
density is evaluated, exactly as a time solver would see them, but with no threading, MPI, file
output or logging around it. The interaction system materials get a smaller body by default.

The step size is part of the fixture because the rotated formulations update their stress
history incrementally, so an arbitrary value would change both the state and the code path.
"""
function material_fixture(mat; n=default_n(mat), m=horizon_ratio(mat), kwargs...)
    body = cube(mat; n, m, kwargs...)
    solver = VelocityVerlet(steps=1)
    dh = handler(body, solver)
    return deform!((chunk=dh.chunks[1], t=solver.Δt, Δt=solver.Δt))
end

default_n(::Peridynamics.AbstractMaterial) = 6
default_n(::Peridynamics.AbstractInteractionSystemMaterial) = 5

"""
    deform!(fixture; F, amplitude)

Impose a small non-identity `F` plus an inhomogeneous perturbation on the chunk of `fixture` and
return it.

Both parts matter: an undeformed body has zero bond strain and short-circuits, and a purely
homogeneous deformation makes the zero-energy mode stabilization of the correspondence
formulations vanish, which is the expensive part of exactly those materials. The displacement is
treated as an increment over the step size, so the velocity fields are populated too, which
`CRMaterial` and `RKCRMaterial` need for the velocity gradient.
"""
function deform!(fixture; F=[1.01 0.003 0.0; 0.0 0.995 0.0; 0.0 0.0 1.0], amplitude=1e-4)
    (; chunk, Δt) = fixture
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
            velocity[d, i] = u / Δt
            velocity_half[d, i] = u / Δt
        end
    end
    return fixture
end

"""
    apply_deformation!(chunk, F, Δt)

Set the current positions of all points of `chunk` to `F * X` and the velocities to the
corresponding increment over `Δt`, i.e. a homogeneous deformation with deformation gradient `F`.
"""
function apply_deformation!(chunk, F, Δt)
    (; system, storage) = chunk
    for i in Peridynamics.each_point_idx(system)
        Xi = Peridynamics.get_vector(system.position, i)
        xi = F * Xi
        Peridynamics.update_vector!(storage.position, i, xi)
        vi = (xi - Xi) / Δt
        Peridynamics.update_vector!(storage.velocity, i, vi)
        Peridynamics.update_vector!(storage.velocity_half, i, vi)
    end
    return chunk
end

"""
    force_density!(fixture)
    force_density!(chunk, t, Δt)

Run the complete force density calculation, the hot loop of every simulation.

The reproducing kernel materials override `calc_force_density!` for the data handler and not for
the chunk, because they need a halo exchange between the gradient weights and the forces.
Calling the chunk method alone would silently skip `calc_weights_and_defgrad!` and check half
the work, so those materials get their own method here.
"""
force_density!(fixture) = force_density!(fixture.chunk, fixture.t, fixture.Δt)

force_density!(chunk, t, Δt) = Peridynamics.calc_force_density!(chunk, t, Δt)

function force_density!(chunk::Peridynamics.BodyChunk{<:Peridynamics.AbstractBondSystem,
                                                      <:Peridynamics.AbstractRKCMaterial},
                        t, Δt)
    Peridynamics.calc_weights_and_defgrad!(chunk, t, Δt)
    Peridynamics.calc_force_density!(chunk, t, Δt)
    return nothing
end

"""
    gradient_weights!(fixture)

Recompute the gradient weights of a reproducing kernel material for every point of `fixture`.

[`force_density!`](@ref) does not reach this. The weights are only recomputed where damage has
just grown, and `calc_damage!` rewrites the `update_gradients` flag at the start of every force
calculation, so an undamaged body never enters it and the cost stays invisible.
"""
gradient_weights!(fixture) = Peridynamics.initialize!(fixture.chunk)

"Whether [`gradient_weights!`](@ref) does anything for `mat`."
has_gradient_weights(mat) = mat isa Peridynamics.AbstractRKCMaterial

#===== running things =====#

"""
    run_threads(job; n_chunks=1)

`Peridynamics.submit_threads(job, n_chunks)` without any terminal output. `submit` is the only
function that sets the global quiet flag, so a direct `submit_threads` call prints the banner
unless something set the flag before — and which item runs first is an accident of file order.
"""
function run_threads(job; n_chunks=1)
    Peridynamics.set_quiet!(true)
    return Peridynamics.submit_threads(job, n_chunks)
end

#===== condition functions =====#
# Named functions for boundary and initial conditions. `velocity_bc!` and friends specialize on
# the function type, and with coverage on every distinct anonymous function costs a fresh
# compilation of the whole condition machinery. Reusing these keeps it to a handful.

"`f(t)` boundary condition: `1`."
f_one(t) = 1.0
"`f(t)` boundary condition: `t`."
f_t(t) = t
"`f(t)` boundary condition: `2t`."
f_2t(t) = 2t
"`f(p, t)` boundary condition: `p[1] * t`."
f_pt(p, t) = p[1] * t
"`f(p, t)` boundary condition: `p[2] * t`."
f_pt2(p, t) = p[2] * t
"`f(p)` condition: `2`."
f_p(p) = 2.0
"`f(p)` condition: `p[1]`."
f_p1(p) = p[1]
"Invalid: two arguments with the wrong names."
f_bad_ab(a, b) = a * b
"Invalid: three arguments."
f_bad_ktu(k, t, u) = k * t * u
"Invalid: one argument with the wrong name."
f_bad_k(k) = 3.456

#===== randomness =====#

"""
    rng()

A seeded random number generator for the few tests that need random data. Never assert an exact
value that was derived from it: the stream of `Xoshiro` is not guaranteed to be identical
across Julia versions. Prefer the deterministic bodies above.
"""
rng() = Xoshiro(0x5eed)

#===== MPI subprocesses =====#

"""
    run_mpi_script(script, args...; nranks=2, expect_success=true)

Run `test/mpi/scripts/<script>` with `mpiexec -n nranks` in a fresh Julia process that uses the
same project and the same flags as the test process (`Base.julia_cmd()` forwards
`--code-coverage`, so the ranks write coverage files too). Returns whether the exit status is
the expected one; the combined output of all ranks is logged as an error if it is not, so a
failing check inside the script is readable in the test log.
"""
function run_mpi_script(script::AbstractString, args...; nranks::Int=2, expect_success=true)
    path = joinpath(@__DIR__, "..", "mpi", "scripts", script)
    mpiexec = Peridynamics.MPI.mpiexec()
    jlcmd = Base.julia_cmd()
    pdir = pkgdir(Peridynamics)
    script_args = String[string(a) for a in args]
    cmd = `$(mpiexec) -n $(nranks) $(jlcmd) --project=$(pdir) $(path) $(script_args)`
    log = IOBuffer()
    ok = success(pipeline(ignorestatus(cmd); stdout=log, stderr=log))
    if ok != expect_success
        @error "unexpected exit status of MPI script" script expect_success ok output=String(take!(log))
    end
    return ok == expect_success
end
