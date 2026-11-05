# A simulation whose boundary condition throws on one rank. That is not a `NaNError`, which is
# synchronized across the ranks, so `submit` has to abort MPI to avoid a deadlock: the process
# exits non-zero and the error is in the logfile.
#     mpiexec -n 2 julia --project=<Peridynamics> test/mpi/scripts/abort.jl <path>
using Peridynamics
const PATH = ARGS[1]
Δt = 1e-7
pos = [0.0 1.0; 0.0 0.0; 0.0 0.0]
body = Body(BBMaterial(), pos, [1.0, 1.0])
material!(body, horizon=1.5, rho=8000, E=210e9)
velocity_bc!(body, :all_points, :x) do t
    t < 1.9Δt ? 1.0 : error("some weird error occurred!\n")
end
job = Job(body, VelocityVerlet(steps=5, stepsize=Δt); path=PATH)
submit(job; quiet=true) # aborts
