# The benchmark suite. Run it and save a baseline with `benchmark/run.jl`, compare two
# baselines with `benchmark/compare.jl`.
#
# Do not activate an environment here: PkgBenchmark and AirspeedVelocity activate one
# themselves and a `Pkg.activate` in this file breaks them.

using BenchmarkTools
using Peridynamics

include(joinpath(pkgdir(Peridynamics), "jobs", "jobs.jl"))

const SUITE = BenchmarkGroup()

# One entry per material, plus a few single-knob variations of `CMaterial` for the constitutive
# model, the kernel function and the zero-energy mode stabilization. Those three are shared
# code, so varying them on every material would measure the same functions over and over.
#
# A new material only needs a line here to get a force density benchmark. The body it is
# measured on comes from `fixture_size` in `jobs/jobs.jl`.
const MATERIALS = [
    "BB" => BBMaterial(),
    "BB-ESC" => BBMaterial{EnergySurfaceCorrection}(),
    "DHBB" => DHBBMaterial(),
    "GBB" => GBBMaterial(),
    "OSB" => OSBMaterial(),
    "OSB-ESC" => OSBMaterial{EnergySurfaceCorrection}(),
    "C-ZEMSilling" => CMaterial(model=SaintVenantKirchhoff(), zem=ZEMSilling()),
    "C-ZEMWan" => CMaterial(model=SaintVenantKirchhoff(), zem=ZEMWan()),
    "CR" => CRMaterial(),
    "BAC" => BACMaterial(),
    "RKC-C1" => RKCMaterial(monomial=:C1),
    "RKC-RK1" => RKCMaterial(monomial=:RK1),
    "RKCR" => RKCRMaterial(),
    "CKI" => CKIMaterial(),
]

# The hot loop, measured on a single chunk so that threading, MPI, file output and logging stay
# out of the measured region. A fresh sample for every call, and `evals=1` so that a sample
# cannot contain more than one update. The memory column of `compare.jl` is where an allocating
# force density calculation shows up; it should read `0 bytes` for every material.
SUITE["force_density"] = BenchmarkGroup()
for (name, mat) in MATERIALS
    fixture = material_fixture(mat)
    SUITE["force_density"][name] = @benchmarkable(force_density!(sample),
                                                  setup=(sample=deepcopy($fixture)), evals=1)
end

# A separate group, because `force_density!` never reaches this path on an undamaged body and
# it would otherwise be invisible. See `gradient_weights!` in `jobs/jobs.jl`.
SUITE["gradient_weights"] = BenchmarkGroup()
for (name, mat) in MATERIALS
    has_gradient_weights(mat) || continue
    fixture = material_fixture(mat)
    SUITE["gradient_weights"][name] = @benchmarkable gradient_weights!($fixture)
end

# Whole simulations. None of these jobs is given a `path`, so nothing is written to disk inside
# the measured region, and `quiet=true` keeps the logging out of it too.
SUITE["job"] = BenchmarkGroup()
for (name, mat) in ("BB" => BBMaterial(), "OSB" => OSBMaterial(), "C" => CMaterial())
    SUITE["job"]["mode_i $name"] = @benchmarkable(submit(mode_i(; mat=$mat, npxy=20,
                                                                time=6e-5); quiet=true),
                                                  seconds=20)
    SUITE["job"]["wave_in_bar $name"] = @benchmarkable(submit(wave_in_bar(; mat=$mat,
                                                                         time=4e-6);
                                                              quiet=true), seconds=20)
end
