# MPI scaling benchmark: how the wall time of a single fixed problem falls as ranks are added,
# and where it stops falling.
#
#     mpiexecjl -n <ranks> julia --project=benchmark/hpc benchmark/hpc/scaling.jl [options]
#
# Options, all with `--key=value`:
#
#     --size=smoke|small|large   problem size, see `SIZES` below  (default: smoke)
#     --mat=BB                   material, see `MATS` below       (default: BB)
#     --label=<name>             subdirectory of `results/`       (default: the host name)
#     --out=<path>               where `results/` lives           (default: next to this)
#     --force                    run `large` without MPI anyway
#
# One invocation measures one point of the scaling curve and writes one CSV file. Run it once
# per rank count, then draw the curve with `benchmark/hpc/report.jl`.
#
# This file does no package management: under `mpiexec` every rank would resolve and precompile
# the same environment at the same moment, which is a good way to corrupt a depot on a shared
# file system. `submit.sh` instantiates once, serially, before the launcher is called.

using Peridynamics
using Printf

include(joinpath(pkgdir(Peridynamics), "jobs", "jobs.jl"))

# The measured problem is `wave_in_bar`: a slender bar is the shape that makes a
# one-dimensional decomposition meaningful, and it is the problem of the reference scaling
# study. `m = 2.015` is from that study as well, and smaller than the 3.015 used elsewhere
# here, so more of the time goes into the bonds and less into the neighbor bookkeeping.
const HORIZON_RATIO = 2.015

# `npyz` is the number of points across the bar, so the point count grows with its cube.
# `large` is 2.7 million points, the size of the reference study.
const SIZES = (smoke=(npyz=6, steps=50),
               small=(npyz=15, steps=200),
               large=(npyz=30, steps=1000))

const MATS = Dict("BB" => BBMaterial(), "OSB" => OSBMaterial(), "C" => CMaterial(),
                  "CR" => CRMaterial(), "RKC" => RKCMaterial())

const RESULT_HEADER = "ranks,threads,seconds,points,steps,material,npyz,m,julia,commit,host,date"

function usage()
    msg = "usage: mpiexecjl -n <ranks> julia --project=benchmark/hpc " *
          "benchmark/hpc/scaling.jl [options]\n"
    msg *= "  --size=$(join(keys(SIZES), "|"))\n"
    msg *= "  --mat=$(join(sort(collect(keys(MATS))), "|"))\n"
    msg *= "  --label=<name>   --out=<path>   --force\n"
    return msg
end

function parse_options(args)
    opts = Dict("size" => "smoke", "mat" => "BB", "label" => hostname(),
                "out" => joinpath(@__DIR__, "results"), "force" => "false")
    for arg in args
        m = match(r"^--([a-z]+)(?:=(.*))?$", arg)
        isnothing(m) && error("cannot parse `$(arg)`.\n" * usage())
        haskey(opts, m[1]) || error("unknown option `--$(m[1])`.\n" * usage())
        opts[m[1]] = isnothing(m[2]) ? "true" : m[2]
    end
    haskey(SIZES, Symbol(opts["size"])) || error("unknown size `$(opts["size"])`.\n" * usage())
    haskey(MATS, opts["mat"]) || error("unknown material `$(opts["mat"])`.\n" * usage())
    return opts
end

hostname() = something(get(ENV, "SLURM_CLUSTER_NAME", nothing), gethostname())

"""
    launched_by_mpi()

Whether this process was started by an MPI launcher, judged from the environment it set.

Peridynamics falls back to its threads path when it sees a single rank. That is the right
default and the wrong thing here: a single rank is a legitimate point of a scaling curve, and
measuring it on a different code path than the rest of the curve makes the speedups meaningless.
So the launcher is detected directly and the MPI path is forced.
"""
function launched_by_mpi()
    return any(k -> startswith(k, "PMI_") || startswith(k, "OMPI_") ||
                    k == "MPI_LOCALRANKID" || k == "SLURM_PROCID", keys(ENV))
end

"The commit of the Peridynamics checkout being measured, or `unknown` outside a git tree."
function commit()
    try
        return readchomp(`git -C $(pkgdir(Peridynamics)) rev-parse --short HEAD`)
    catch
        return "unknown"
    end
end

"""
    scaling_job(mat, problem)

The job whose wall time is the measurement.

No `path` is given, so nothing is exported inside the measured region: on a parallel file system
the write would dominate the thing being measured. `sort_points` is essential here and the
reason it exists — see [`wave_in_bar`](@ref).
"""
function scaling_job(mat, problem)
    return wave_in_bar(; mat, npyz=problem.npyz, m=HORIZON_RATIO, steps=problem.steps,
                       sort_points=true)
end

"""
    measure(mat, problem)

Wall time in seconds of one run of [`scaling_job`](@ref), measured on the root rank.

The ranks are synchronized before the clock starts, and they have exchanged halos `steps` times
by the time the last one leaves `submit`, so the root rank's time is the time of the slowest
rank to within one exchange.
"""
function measure(mat, problem)
    job = scaling_job(mat, problem)
    mpi_barrier()
    seconds = @elapsed submit(job; quiet=true)
    return seconds, job.spatial_setup.n_points
end

"""
    warmup(mat, problem)

One short run whose time is thrown away, so that the measured run does not carry the
compilation of the force density calculation, the time solver and the halo exchange. That is
easily several seconds — the same order as the measurement itself at high rank counts, and it
would fall unevenly across the configurations being compared.

It uses the same body as the measurement, because a smaller one would decompose into chunks of
a few points each at high rank counts.
"""
function warmup(mat, problem)
    @mpiroot @info "warmup run"
    measure(mat, (npyz=problem.npyz, steps=10))
    return nothing
end

function write_result(opts, ranks, threads, seconds, points, problem)
    dir = joinpath(opts["out"], opts["label"])
    mkpath(dir)
    file = joinpath(dir, @sprintf("%s_%s_n%05d-t%03d.csv", opts["mat"], opts["size"], ranks,
                                  threads))
    isfile(file) && @warn "overwriting an earlier result" file
    open(file, "w") do io
        println(io, RESULT_HEADER)
        println(io, join((ranks, threads, seconds, points, problem.steps, opts["mat"],
                          problem.npyz, HORIZON_RATIO, VERSION, commit(), gethostname(),
                          now_string()), ","))
    end
    @info "result written" file seconds
    return nothing
end

now_string() = Libc.strftime("%Y-%m-%dT%H:%M:%S", time())

function main(args)
    opts = parse_options(args)
    problem = SIZES[Symbol(opts["size"])]
    mat = MATS[opts["mat"]]
    launched_by_mpi() && force_mpi_run!()
    ranks = Peridynamics.mpi_run() ? Peridynamics.mpi_nranks() : 1
    threads = Threads.nthreads()

    # `smoke` and `small` stay available everywhere, because being able to test the script here
    # is what keeps a cluster job from failing after it has waited in a queue.
    if opts["size"] == "large" && !Peridynamics.mpi_run() && opts["force"] == "false"
        error("the `large` size is meant for a cluster and is not being launched under MPI.\n" *
              "Use `--size=smoke` to check the script here, or pass `--force` if this really " *
              "is the machine you want to measure.\n")
    end

    @mpiroot @info "MPI scaling benchmark" size=opts["size"] material=opts["mat"] ranks threads

    warmup(mat, problem)
    @mpiroot @info "measured run"
    seconds, points = measure(mat, problem)
    @mpiroot write_result(opts, ranks, threads, seconds, points, problem)
    return nothing
end

main(ARGS)
