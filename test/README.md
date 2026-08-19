# Test suite

Test items (`@testitem`, TestItemRunner) organized by what they check. Run them with

```
julia -t 6 test/runtestitems.jl                       # everything (incl. :skipci)
julia -t 6 test/runtestitems.jl file:test/unit/core/test_study.jl
julia -t 6 test/runtestitems.jl tag:simulation        # every time loop in the suite
julia -t 8 test/runtestitems.jl tag:verification      # ~1 h, never on CI
julia -t 2 -e 'using Pkg; Pkg.test()'                 # what CI runs (PERIDYNAMICS_TESTS=default)
```

## Layout

| Directory | Contents | Tags | CI job |
|---|---|---|---|
| `unit/<dir>/` | unit tests mirroring `src/`: `unit/<dir>/test_<file>.jl` tests `src/<dir>/<file>.jl` | — | `test` |
| `materials/` | one item per property, looping `Fixtures.MATERIALS` (storage fields, force density invariants, stress/energy exports, the material interface) | — / `:slow` | `test` / `extras` |
| `simulation/` | the per-commit time loops: symmetry of every solver, uniform tension | `:simulation` | `test` |
| `mpi/` | three `mpiexec` runs (core paths, abort, threads comparison), scripts in `mpi/scripts/` | `:mpi` | `extras` |
| `perf/` | allocation and type stability checks | `:perf` | `extras` |
| `quality/` | Aqua and the convention guards | `:lint` | `extras` |
| `verification/` | physics against closed forms and convergence rates | `:verification, :skipci` | never |
| `setup/` | the shared `@testmodule`s: `Fixtures` (materials table, bodies, chunks, condition functions, `rng()`) and `TestMaterialImpl` | | |

`test/runtests.jl` selects by `PERIDYNAMICS_TESTS` (`default`, `extras`, `all`); see
`.github/workflows/CI.yml`.

## Conventions

Enforced by `quality/test_conventions.jl` where possible:

- one item, one behaviour; name `"<function or type>: <behaviour>"`;
- every item that runs a time solver carries `:simulation`, every `mpiexec` spawn `:mpi`;
  the per-commit simulations are an explicit allowlist in the guard;
- no per-material copy-paste: loop `Fixtures.MATERIALS` / `MODEL_MATERIALS` in a `@testset`;
- assert derived values; a literal constant needs a comment saying where it comes from;
- no `include`, no printing, no unseeded `rand` inside items (use `Fixtures.rng()`);
  expected warnings are asserted with `@test_logs`; `submit` is called with `quiet=true`;
- `:simulation` items are sized for the code path, accuracy claims live in `verification/`.
