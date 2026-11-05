# Benchmarks

Two separate things live here, and they answer different questions.

| | question | where it runs |
|---|---|---|
| `run.jl` + `compare.jl` | did this change make the code slower? | workstation |
| `hpc/` | how far does one simulation scale across ranks? | cluster |

The problems themselves are defined in `jobs/` at the root of the repository, so that the
benchmarks and the verification tests measure the same simulation.

## Comparing two revisions

```
julia --project=benchmark benchmark/run.jl
julia --project=benchmark benchmark/compare.jl <old> <new>
```

`run.jl` names the baseline after the current commit, with `-dirty` appended when the working
tree is modified, and writes it to `benchmark/baselines/` (not tracked by git). Pass a name as
an argument to override.

`compare.jl` prints one section per group, slowest benchmark first, and dims everything inside
the 10 % tolerance. The marker at the end of each row repeats the verdict as a character — `▲`
slower, `▼` faster, `·` within the tolerance.

The memory column shows the amount allocated, and switches to a relative change once that amount
moves. The `force_density` rows should read `0 bytes`, which `test/perf/perf.jl` asserts as
well: a force density calculation that starts allocating is a defect and not a trade-off to be
weighed against a speedup.

The groups are `force_density` (the hot loop, `calc_force_density!` on one chunk),
`gradient_weights` (the reproducing kernel materials only) and `job` (whole simulations).

Because the suite follows the `const SUITE` convention it also works with AirspeedVelocity:

```
benchpkg Peridynamics --rev=main,HEAD
```

## MPI scaling on a cluster

`hpc/` measures one fixed problem under an increasing number of MPI ranks. Nothing in it is
reachable from `run.jl`, and the large size refuses to start outside MPI, so it cannot cost
anything on a workstation by accident.

```
sbatch benchmark/hpc/submit.sh            # adjust the SBATCH header once per cluster
```

The submit script instantiates the environment once and serially, then runs one invocation per
point of the scaling curve. Each invocation writes one CSV file into `hpc/results/<machine>/`.
Copy that directory back and draw the curve with

```
julia --project=benchmark/hpc benchmark/hpc/report.jl
```

Before queueing anything, check the script itself locally. This takes seconds:

```
mpiexecjl -n 4 julia --project=benchmark/hpc benchmark/hpc/scaling.jl --size=smoke
```
