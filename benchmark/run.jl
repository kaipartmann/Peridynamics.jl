# Run the benchmark suite and save the result as a baseline.
#
#     julia --project=benchmark benchmark/run.jl [name]
#
# The baseline is named after the current git commit unless a name is given, so the usual
# workflow is to check out `main`, run this, check out the branch, run it again, and compare
# the two with `benchmark/compare.jl`.

import Pkg

# `Manifest.toml` is not committed, so bind Peridynamics to this checkout explicitly: otherwise
# switching revisions would silently keep benchmarking a registered release.
Pkg.activate(@__DIR__)
Pkg.develop(Pkg.PackageSpec(path=normpath(joinpath(@__DIR__, ".."))))

using BenchmarkTools

include(joinpath(@__DIR__, "benchmarks.jl"))

const BASELINES = joinpath(@__DIR__, "baselines")

function default_name()
    sha = readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
    dirty = !isempty(readchomp(`git -C $(@__DIR__) status --porcelain`))
    return dirty ? sha * "-dirty" : sha
end

function main(args)
    name = isempty(args) ? default_name() : args[1]
    mkpath(BASELINES)
    file = joinpath(BASELINES, name * ".json")
    isfile(file) && @warn "overwriting existing baseline" file
    results = run(SUITE; verbose=true)
    BenchmarkTools.save(file, results)
    @info "baseline saved" file
    return nothing
end

main(ARGS)
