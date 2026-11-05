# Job creators for the `Study` tests.

"""
    tiny_job(setup, root)

The job creator every `Study` test uses: a ten-point line body with `E = setup.E`, a velocity
initial condition and a `VelocityVerlet` run of `setup.n_steps` steps, exported every 5 steps
to `joinpath(root, "sim_<E>_<n_steps>")`. Sharing one creator keeps `Study{typeof(creator)}`
and everything it calls compiled once instead of once per item.

`setup` may carry a `fail::Bool` field; a failing job gets a boundary condition that throws in
the first time step, so the job is created fine and fails in `submit` — an error that does not
depend on the operating system or on file permissions.
"""
function tiny_job(setup::NamedTuple, root::String)
    body = Body(BBMaterial(), line10_position(), ones(10))
    material!(body; horizon=1.5, E=setup.E, rho=1, Gc=1)
    velocity_ic!(body, :all_points, :x, 1.0)
    if get(setup, :fail, false)
        velocity_bc!(failing_bc, body, :all_points, :y)
    end
    vv = VelocityVerlet(steps=setup.n_steps)
    path = joinpath(root, "sim_$(setup.E)_$(setup.n_steps)")
    return Job(body, vv; path=path, freq=5)
end

"A boundary condition function that throws, to make a job fail deliberately."
failing_bc(t) = error("deliberately failing job")

function line10_position()
    position = zeros(3, 10)
    position[1, :] .= 0.0:9.0
    return position
end

"Setups for `tiny_job`: three jobs with different stiffness and step counts."
const TINY_SETUPS = [(; E=1.0, n_steps=5), (; E=2.0, n_steps=10), (; E=3.0, n_steps=15)]
