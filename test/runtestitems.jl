#!/usr/bin/env julia
#
# Run Peridynamics.jl test items from the command line.
#
# USAGE:
#   julia -t 6 test/runtestitems.jl [--list] [selector ...]
#
# SELECTORS (a test item runs only if every selector matches; no selector = run all):
#   <name fragment>   case-insensitive substring of the test-item name (AND-ed)
#   file:<path>       test item is defined in this file (OR-ed across file selectors)
#   tag:<tag>         test item carries this tag (AND-ed)
#
# TAGS used in this suite (see test/runtests.jl for what CI runs):
#   simulation    the item runs a time solver; `tag:simulation` lists every simulation
#   mpi           the item spawns `mpiexec`; runs in the `mpi` CI job
#   perf          allocation and type-stability checks; `extras` CI job
#   lint          Aqua and the suite convention guards; `extras` CI job
#   slow          heavy but valuable items moved to the `extras` CI job
#   verification  physics checked against closed forms or convergence rates; never in CI
#   skipci        never runs in CI
#
# EXAMPLES:
#   Run all:       julia -t 6 test/runtestitems.jl
#   By name:       julia -t 6 test/runtestitems.jl BBMaterial dynamic
#   By file:       julia -t 6 test/runtestitems.jl file:test/unit/core/test_halo_exchange.jl
#   By tag:        julia -t 6 test/runtestitems.jl tag:mpi
#   Simulations:   julia -t 6 test/runtestitems.jl tag:simulation
#   Verification:  julia -t 8 test/runtestitems.jl tag:verification        # ~1 h
#                  julia -t 8 test/runtestitems.jl tag:verification wave   # only the wave items
#   Combined:      julia -t 6 test/runtestitems.jl BBMaterial file:test/unit/core/test_halo_exchange.jl tag:mpi
#   Preview:       julia -t 6 test/runtestitems.jl --list tag:mpi
#
# Selectors that match nothing are an error, so a typo cannot look like a green run.
#
# REQUIREMENTS: TestEnv.jl in the global environment, install it once with
#   julia -e 'import Pkg; Pkg.add("TestEnv")'
#
# SANDBOXED ENVIRONMENTS (e.g. containers):
#   The `tag:mpi` test items spawn `mpiexec -n 2 julia ...` subprocesses, which need a
#   writable Julia depot and permission to start local processes. If the default depot
#   (~/.julia) is read-only, this script prepends a writable scratch depot and keeps the
#   read-only depot as a read fallback, so packages and artifacts (including the MPI
#   binaries) are still found and nothing has to be downloaded. JULIA_DEPOT_PATH is
#   exported so the mpiexec subprocesses inherit the same setup.
#   Set PERIDYNAMICS_TEST_DEPOT=<path> to force a specific scratch depot. It must live
#   outside this repository, otherwise TestItemRunner would scan it for test items.
#   If MPI still fails, the sandbox blocks process spawning or local sockets rather than
#   file access; run everything except MPI with `--list`-verified selectors, e.g.
#   `julia -t 6 test/runtestitems.jl` after checking which items carry `tag:mpi`.

import Pkg

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEPOT_SEP = Sys.iswindows() ? ';' : ':'

"""
    is_writable(dir)

Check whether `dir` exists (or can be created) and accepts writes, without leaving
anything behind. Used to detect read-only depots in sandboxed environments.
"""
function is_writable(dir::AbstractString)
    probe = joinpath(dir, ".runtestitems-write-probe-$(getpid())")
    try
        mkpath(dir)
        touch(probe)
        rm(probe; force=true)
        return true
    catch
        return false
    end
end

"""
    setup_depot!()

Make sure `DEPOT_PATH[1]` is writable, because that is where Julia stores precompilation
caches. If the default depot is read-only, prepend a writable scratch depot and keep the
read-only depot as a read fallback for already installed packages and artifacts.
"""
function setup_depot!()
    forced = get(ENV, "PERIDYNAMICS_TEST_DEPOT", "")
    if isempty(forced) && !isempty(DEPOT_PATH) && is_writable(first(DEPOT_PATH))
        return nothing # default depot is fine, nothing to do
    end
    candidates = isempty(forced) ? [joinpath(tempdir(), "peridynamics-test-depot")] :
                 [forced]
    for candidate in candidates
        if startswith(normpath(abspath(candidate)), REPO_ROOT)
            error("The scratch depot must be outside the repository, but it is at:\n" *
                  "    $(candidate)\n" *
                  "TestItemRunner scans the whole repository for test items and would " *
                  "descend into the depot.")
        end
    end
    scratch = findfirst(is_writable, candidates)
    if scratch === nothing
        default = isempty(DEPOT_PATH) ? "<none>" : first(DEPOT_PATH)
        error("No writable Julia depot available. The default depot is\n" *
              "    $(default)\n" *
              "and none of these fallbacks could be written to:\n" *
              join("    " .* candidates, "\n") * "\n" *
              "Set PERIDYNAMICS_TEST_DEPOT=<writable path outside the repository>.")
    end
    depot = candidates[scratch]
    depot in DEPOT_PATH || pushfirst!(DEPOT_PATH, depot)
    # exported so the mpiexec subprocesses of the `tag:mpi` test items inherit it
    ENV["JULIA_DEPOT_PATH"] = join(DEPOT_PATH, DEPOT_SEP)
    reason = isempty(forced) ? "the default depot is read-only" :
             "PERIDYNAMICS_TEST_DEPOT is set"
    @info "Writing to a scratch depot, because $(reason)" depot
    return nothing
end

"""
    parse_selectors(args)

Split the command line into name, file and tag selectors, plus the `--list` flag.
"""
function parse_selectors(args)
    list_only = false
    names = String[]
    files = String[]
    tags = String[]
    for argument in args
        if argument in ("--list", "-l")
            list_only = true
        elseif startswith(lowercase(argument), "tag:")
            tag = argument[5:end]
            isempty(tag) && error("A tag selector must use the form tag:<tag>.")
            push!(tags, lowercase(tag))
        elseif startswith(lowercase(argument), "file:")
            file = argument[6:end]
            isempty(file) && error("A file selector must use the form file:<path>.")
            push!(files, resolve_test_file(file))
        elseif startswith(argument, "-")
            error("Unknown option `$(argument)`. Supported options: --list.")
        else
            push!(names, lowercase(argument))
        end
    end
    return (; list_only, names, files, tags)
end

"""
    resolve_test_file(file)

Turn a file selector into an absolute path, accepting paths relative to the current
directory as well as relative to the repository root.
"""
function resolve_test_file(file::AbstractString)
    path = normpath(abspath(file))
    isfile(path) && return path
    alternative = normpath(joinpath(REPO_ROOT, file))
    isfile(alternative) && return alternative
    error("The file selector `file:$(file)` does not point to an existing file. Tried:\n" *
          "    $(path)\n    $(alternative)")
end

"""
    matches(ti, selectors)

Check a test item (a named tuple with `name`, `filename` and `tags`) against all
selectors. Name and tag selectors are AND-ed, file selectors are OR-ed.
"""
function matches(ti, selectors)
    name = lowercase(ti.name)
    all(contains(name, filter) for filter in selectors.names) || return false
    if !isempty(selectors.files)
        normpath(abspath(ti.filename)) in selectors.files || return false
    end
    tags = [lowercase(String(tag)) for tag in ti.tags]
    all(in(tags), selectors.tags) || return false
    return true
end

setup_depot!()

Pkg.activate(REPO_ROOT; io=devnull)

if Base.find_package("TestEnv") === nothing
    error("TestEnv.jl is not installed in Julia's global environment. Install it with:" *
          "\n\n    julia -e 'import Pkg; Pkg.add(\"TestEnv\")'")
end

using TestEnv

TestEnv.activate()

using TestItemRunner

const SELECTORS = parse_selectors(ARGS)
const N_MATCHED = Ref(0)

"""
    select(ti)

Filter callback for TestItemRunner, called exactly once per discovered test item. It
counts the matches and, in `--list` mode, prints them instead of running them.
"""
function select(ti)
    matches(ti, SELECTORS) || return false
    N_MATCHED[] += 1
    if SELECTORS.list_only
        tags = isempty(ti.tags) ? "" : "  tags: " * join(ti.tags, ", ")
        println(rpad(ti.name, 60), "  ", relpath(ti.filename, REPO_ROOT), tags)
        return false
    end
    return true
end

@run_package_tests verbose=true filter=select

if N_MATCHED[] == 0
    error("No test item matched the selectors: $(join(ARGS, " "))")
elseif SELECTORS.list_only
    println("\n$(N_MATCHED[]) test item(s) matched.")
end
