# Turn the CSV files written by `scaling.jl` into a scaling table.
#
#     julia --project=benchmark/hpc benchmark/hpc/report.jl [label] [--out=<path>]
#
# `label` is the subdirectory of `results/` to report on, i.e. the machine. With no argument
# every label found is reported. This only reads files and prints, so it runs anywhere,
# including on the workstation after the results have been copied back from the cluster.

using Printf

include(joinpath(@__DIR__, "..", "pretty.jl"))

const RESULTS = Ref(joinpath(@__DIR__, "results"))

const W_RANKS, W_THREADS, W_CORES, W_TIME, W_SPEEDUP, W_EFF, W_BAR = 7, 8, 7, 12, 9, 12, 12

function read_result(file)
    lines = readlines(file)
    length(lines) == 2 || error("`$(file)` is not a result file written by scaling.jl")
    header = split(lines[1], ',')
    values = split(lines[2], ',')
    r = Dict(zip(header, values))
    return (ranks=parse(Int, r["ranks"]), threads=parse(Int, r["threads"]),
            seconds=parse(Float64, r["seconds"]), points=parse(Int, r["points"]),
            steps=parse(Int, r["steps"]), material=r["material"], npyz=parse(Int, r["npyz"]),
            commit=r["commit"])
end

function read_label(label)
    dir = joinpath(RESULTS[], label)
    isdir(dir) || error(missing_label_message(label, dir))
    files = sort(filter(f -> endswith(f, ".csv"), readdir(dir; join=true)))
    isempty(files) && error("no result files in $(dir)")
    return map(read_result, files)
end

function missing_label_message(label, dir)
    msg = "no results labelled `$(label)`!\n"
    msg *= "Looked in $(dir).\n"
    available = isdir(RESULTS[]) ? sort(filter(f -> isdir(joinpath(RESULTS[], f)),
                                               readdir(RESULTS[]))) : String[]
    if isempty(available)
        msg *= "There are no results yet. Measure one with\n"
        msg *= "    mpiexecjl -n <ranks> julia --project=benchmark/hpc " *
               "benchmark/hpc/scaling.jl --size=smoke\n"
        msg *= "or copy the result files back from the cluster.\n"
    else
        msg *= "Available labels:\n"
        for a in available
            msg *= "    $(a)\n"
        end
    end
    return msg
end

"""
    scaling_group(r)

The set of results `r` may be compared against. A speedup only means something within one
problem, i.e. the same material at the same size.
"""
function scaling_group(r)
    return string(r.material, "   ", prettycount(r.points), " points, ",
                  prettycount(r.steps), " steps")
end

"""
    efficiency_style(eff)

Color for a parallel efficiency. The thresholds are a reading aid and not a verdict on the code:
below about half, the added cores are mostly waiting for each other.
"""
function efficiency_style(eff)
    eff ≥ 0.75 && return (color=:green, bold=false)
    eff ≥ 0.50 && return (color=:yellow, bold=false)
    return (color=:red, bold=true)
end

"A bar of `W_BAR` cells, filled in proportion to `eff`."
bar(eff) = rpad("█"^clamp(round(Int, W_BAR * eff), 0, W_BAR), W_BAR)

prettyseconds(s) = s < 60 ? @sprintf("%.2f s", s) : @sprintf("%d:%04.1f min", s ÷ 60, s % 60)

function print_group(results)
    sort!(results; by=r -> (r.ranks * r.threads, r.ranks))
    ref = first(results)
    colheads("  ", lpad("ranks", W_RANKS), lpad("threads", W_THREADS), lpad("cores", W_CORES),
             lpad("wall time", W_TIME), lpad("speedup", W_SPEEDUP),
             lpad("efficiency", W_EFF), "  scaling")
    for r in results
        cores = r.ranks * r.threads
        speedup = ref.seconds / r.seconds
        # against the smallest configuration present and not against one core, because a problem
        # of this size usually does not fit on one core at all
        eff = speedup / (cores / (ref.ranks * ref.threads))
        style = efficiency_style(eff)
        @printf("  %*d%*d%*d", W_RANKS, r.ranks, W_THREADS, r.threads, W_CORES, cores)
        print(lpad(prettyseconds(r.seconds), W_TIME))
        print(lpad(@sprintf("%.2f ×", speedup), W_SPEEDUP))
        printstyled(lpad(@sprintf("%.0f %%", 100 * eff), W_EFF); style.color, style.bold)
        printstyled("  ", bar(eff), "\n"; style.color)
    end
    return nothing
end

function report(label)
    results = read_label(label)
    commits = unique(r.commit for r in results)
    subtitle = @sprintf("%d run%s · commit %s", length(results),
                        isone(length(results)) ? "" : "s", join(commits, ", "))
    title("MPI scaling   $(label)", subtitle)
    # a table mixing two commits is not a scaling study, it is two of them overlaid
    length(commits) > 1 &&
        printstyled(" these results come from different commits!\n"; color=:yellow, bold=true)
    for group in sort(unique(scaling_group.(results)))
        section(group)
        print_group(filter(r -> scaling_group(r) == group, results))
    end
    println()
    return nothing
end

function main(args)
    out = findfirst(a -> startswith(a, "--out="), args)
    isnothing(out) || (RESULTS[] = args[out][7:end])
    labels = filter(a -> !startswith(a, "--"), args)
    if isempty(labels)
        isdir(RESULTS[]) || error(missing_label_message("", RESULTS[]))
        labels = sort(filter(f -> isdir(joinpath(RESULTS[], f)), readdir(RESULTS[])))
        isempty(labels) && error(missing_label_message("", RESULTS[]))
    end
    foreach(report, labels)
    return nothing
end

main(ARGS)
