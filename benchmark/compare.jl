# Compare two baselines saved by `benchmark/run.jl`.
#
#     julia --project=benchmark benchmark/compare.jl <old> <new>
#
# Names are the baseline names without the `.json`. The report is grouped the way the suite is
# grouped, and the slowest benchmark of a group comes first.
#
# The tolerance is 10 percent and not the 5 percent BenchmarkTools defaults to, because below
# that a difference on a normal machine is noise.
import Pkg
Pkg.activate(@__DIR__)

using BenchmarkTools
using Printf

include(joinpath(@__DIR__, "pretty.jl"))

const TIME_TOLERANCE = 0.10
const BASELINES = joinpath(@__DIR__, "baselines")

# Column widths. Every numeric column is two characters wider than its widest entry so that
# adjacent numbers keep a gap, and `W_MARKER` is wider still so that the verdict marker does not
# read as part of the memory column beside it.
const W_NAME, W_TIME, W_CHANGE, W_MEM, W_MARKER = 18, 13, 12, 12, 3

function load_baseline(name)
    file = joinpath(BASELINES, name * ".json")
    isfile(file) || error(missing_baseline_message(name, file))
    return only(BenchmarkTools.load(file))
end

function missing_baseline_message(name, file)
    msg = "no baseline named `$(name)`!\n"
    msg *= "Looked for $(file).\n"
    available = isdir(BASELINES) ?
                sort([first(splitext(f)) for f in readdir(BASELINES) if endswith(f, ".json")]) :
                String[]
    if isempty(available)
        msg *= "There are no saved baselines yet. Create one with\n"
        msg *= "    julia --project=benchmark benchmark/run.jl [name]\n"
    else
        msg *= "Available baselines:\n"
        for a in available
            msg *= "    $(a)\n"
        end
    end
    return msg
end

"The key paths of every leaf benchmark of `group`, as `\"a / b\"` strings."
leaf_names(group) = Set(join(keys, " / ") for (keys, _) in BenchmarkTools.leaves(group))

"The leaf of `group` at the key path `keys`."
at(group, keys) = foldl(getindex, keys; init=group)

"""
    rows_of(judgement, old_min, new_min)

One named tuple per benchmark, carrying everything the table shows. The group is kept separate
from the name so that the report can be laid out the way the suite is, instead of interleaving
single kernel calls with whole simulations.
"""
function rows_of(judgement, old_min, new_min)
    return map(BenchmarkTools.leaves(judgement)) do (keys, leaf)
        return (group=first(keys), name=join(keys[2:end], " / "),
                t_old=time(at(old_min, keys)), t_new=time(at(new_min, keys)),
                m_old=memory(at(old_min, keys)), m_new=memory(at(new_min, keys)),
                t_ratio=ratio(leaf).time, m_ratio=ratio(leaf).memory,
                verdict=time(leaf), mem_verdict=memory(leaf))
    end
end

"""
    group_order(groups)

The groups of the suite in reading order: the force density first, because it is what the other
two are interpreted against, and whole simulations last, because they are the noisiest. A group
not named here is appended alphabetically.
"""
function group_order(groups)
    preferred = ["force_density", "gradient_weights", "job"]
    known = [g for g in preferred if g in groups]
    return vcat(known, sort([g for g in groups if !(g in known)]))
end

function print_table(rows)
    colheads("  ", rpad("benchmark", W_NAME), lpad("old", W_TIME), lpad("new", W_TIME),
             lpad("change", W_CHANGE), lpad("memory", W_MEM))
    for group in group_order(unique(row.group for row in rows))
        section(group)
        in_group = sort(filter(row -> row.group == group, rows); by=row -> row.t_ratio,
                        rev=true)
        foreach(print_row, in_group)
    end
    return nothing
end

function print_row(row)
    style = verdict_style(row.verdict)
    print("  ", rpad(row.name, W_NAME))
    printstyled(lpad(BenchmarkTools.prettytime(row.t_old), W_TIME); color=:light_black)
    print(lpad(BenchmarkTools.prettytime(row.t_new), W_TIME))
    printstyled(lpad(pct(row.t_ratio), W_CHANGE); style.color, style.bold)
    print_memory(row)
    printstyled(lpad(style.marker, W_MARKER + 1), "\n"; style.color, style.bold)
    return nothing
end

"""
    print_memory(row)

The memory column: the amount allocated while nothing moves, the relative change once it does.

A force density calculation that starts allocating is a defect and not a slowdown to be weighed
against a speedup, so this column prints the amount itself rather than `+0.0 %`, and a row of
the `force_density` group reads `0 bytes`.
"""
function print_memory(row)
    if row.m_old == row.m_new
        text, style = BenchmarkTools.prettymemory(row.m_new), verdict_style(:invariant)
    else
        # a percentage of zero is meaningless
        text = iszero(row.m_old) ? BenchmarkTools.prettymemory(row.m_new) : pct(row.m_ratio)
        style = verdict_style(row.mem_verdict)
    end
    printstyled(lpad(text, W_MEM); style.color, style.bold)
    return nothing
end

"""
    summary_line(rows)

The counts under the headline, which double as the legend for the marker column: each one
carries the marker it stands for.
"""
function summary_line(rows)
    n(v) = count(row -> row.verdict === v, rows)
    tally(v) = "$(n(v)) $(v === :regression ? "slower" :
                          v === :improvement ? "faster" : "unchanged") " *
               verdict_style(v).marker
    benchmarks = string(length(rows), isone(length(rows)) ? " benchmark" : " benchmarks")
    return join((benchmarks, @sprintf("tolerance %d %%", round(Int, 100 * TIME_TOLERANCE)),
                 tally(:regression), tally(:improvement), tally(:invariant)), " · ")
end

function main(args)
    length(args) == 2 || error("usage: julia --project=benchmark benchmark/compare.jl <old> <new>")
    old, new = load_baseline(args[1]), load_baseline(args[2])
    old_min, new_min = minimum(old), minimum(new)

    # `judge` maps over the *intersection* of the two groups, so a benchmark that was added or
    # removed between the two baselines would quietly disappear from the table below while the
    # report still read "invariant". Name them explicitly.
    old_names, new_names = leaf_names(old_min), leaf_names(new_min)
    only_new = sort(collect(setdiff(new_names, old_names)))
    only_old = sort(collect(setdiff(old_names, new_names)))
    if isempty(intersect(old_names, new_names))
        error("`$(args[1])` and `$(args[2])` have no benchmark in common, so there is " *
              "nothing to compare.\nThey were most likely produced by different versions of " *
              "`benchmark/benchmarks.jl`.\n")
    end

    judgement = judge(new_min, old_min; time_tolerance=TIME_TOLERANCE)
    rows = rows_of(judgement, old_min, new_min)
    title("benchmark comparison   $(args[1]) → $(args[2])", summary_line(rows))
    print_table(rows)
    println()

    report_unmatched("only in $(args[2])", only_new)
    report_unmatched("only in $(args[1])", only_old)
    return nothing
end

function report_unmatched(header, names)
    isempty(names) && return nothing
    printstyled(" ", header, " (not compared):\n"; color=:yellow, bold=true)
    for name in names
        printstyled("   ", name, "\n"; color=:yellow)
    end
    println()
    return nothing
end

main(ARGS)
