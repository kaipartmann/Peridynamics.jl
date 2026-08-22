# Guards for the conventions of this test suite. They parse every test file with the same parser
# TestItemRunner uses, so comments, strings and the heredoc programs of the MPI items cannot fool
# them, and they fail with the file, the item and the offending name.
#
# To opt a single item out of a rule, put a comment `# lint-ok: <reason>` on the line of its
# `@testitem` header.
@testmodule Conventions begin
    using TestItemRunner
    using TestItemRunner: JuliaSyntax
    using TestItemRunner.JuliaSyntax: @K_str, kind, children, SyntaxNode
    using TestItemRunner.TestItemDetection: find_test_detail!, our_range

    const TEST_ROOT = normpath(@__DIR__, "..")
    const SRC_ROOT = normpath(@__DIR__, "..", "..", "src")

    # every `.jl` file below `test/`, relative to `test/`
    function test_files()
        files = String[]
        for (root, _, names) in walkdir(TEST_ROOT), name in names
            endswith(name, ".jl") || continue
            push!(files, relpath(joinpath(root, name), TEST_ROOT))
        end
        return sort!(files)
    end

    # one record per `@testitem`: `file`, `name`, `line`, `tags`, the set of identifiers and
    # macro names used in its body, and the `lint-ok` comment of its header line, if any
    struct Item
        file::String
        name::String
        line::Int
        tags::Vector{Symbol}
        names::Set{Symbol}
        lint_ok::String
    end

    function items(file::AbstractString)
        content = read(joinpath(TEST_ROOT, file), String)
        stream = JuliaSyntax.ParseStream(content; version=VERSION)
        JuliaSyntax.parse!(stream; rule=:all)
        tree = JuliaSyntax.build_tree(SyntaxNode, stream)
        detected, setups, errors = [], [], []
        find_test_detail!(tree, detected, setups, errors)
        isempty(errors) || error("test item definition errors in $(file): $(errors)")
        nodes = Dict{UnitRange{Int},SyntaxNode}()
        collect_testitem_nodes!(nodes, tree)
        result = Item[]
        for ti in detected
            node = nodes[ti.range]
            names = Set{Symbol}()
            collect_names!(names, node[end]) # the begin ... end block
            line = TestItemRunner.compute_line_column(content, first(ti.range)).line
            header = split(content, '\n')[line]
            m = match(r"#\s*lint-ok:\s*(.*)$", header)
            lint_ok = m === nothing ? "" : String(strip(m.captures[1]))
            push!(result, Item(file, ti.name, line, Symbol.(ti.option_tags), names, lint_ok))
        end
        return result
    end

    function collect_testitem_nodes!(nodes, node)
        if kind(node) == K"macrocall" && !JuliaSyntax.is_leaf(node) &&
           node[1].val == Symbol("@testitem")
            nodes[our_range(node)] = node
            return nothing
        end
        JuliaSyntax.is_leaf(node) && return nothing
        for child in children(node)
            collect_testitem_nodes!(nodes, child)
        end
        return nothing
    end

    function collect_names!(names, node)
        if JuliaSyntax.is_leaf(node)
            k = kind(node)
            if k == K"Identifier" || k == K"MacroName"
                val = node.val
                val isa Symbol && push!(names, val)
            end
            return nothing
        end
        for child in children(node)
            collect_names!(names, child)
        end
        return nothing
    end

    # all items of the suite
    all_items() = reduce(vcat, (items(f) for f in test_files()); init=Item[])

    # a readable location for failure messages
    location(it::Item) = "$(it.file):$(it.line) \"$(it.name)\""

    # Functions that start a time loop. An item using any of them is a `:simulation`.
    const TIME_LOOP = Set([:submit, :submit!, :submit_threads, :solve!])
    # Names that are forbidden in test items, with the reason.
    const FORBIDDEN = [
        :include => "define a @testmodule instead",
        :println => "tests must not print; use @debug",
        :print => "tests must not print; use @debug",
        :printstyled => "tests must not print; use @debug",
        Symbol("@printf") => "tests must not print; use @debug",
        Symbol("@show") => "tests must not print; use @debug",
        :rand => "use a deterministic fixture or the seeded `rng = Fixtures.rng()`",
        :randn => "use a deterministic fixture or the seeded `rng = Fixtures.rng()`",
    ]

    # the unit test directories, `test/unit/<dir>`, which mirror `src/<dir>` one-to-one
    const UNIT_ROOT = joinpath(TEST_ROOT, "unit")
    unit_dirs() = [d for d in readdir(UNIT_ROOT) if isdir(joinpath(UNIT_ROOT, d))]
end

@testitem "conventions: time loops are tagged" tags=[:lint] setup=[Conventions] begin # lint-ok: guard
    # Every item that runs a time solver carries `:simulation`, so that
    #     julia -t 6 test/runtestitems.jl tag:simulation
    # is the complete inventory of simulations in the suite. Items spawning `mpiexec` are `:mpi`.
    untagged_loops = String[]
    untagged_mpi = String[]
    for it in Conventions.all_items()
        isempty(it.lint_ok) || continue
        runs_loop = !isempty(intersect(it.names, Conventions.TIME_LOOP))
        spawns_mpi = :mpiexec in it.names
        tagged = :simulation in it.tags || :mpi in it.tags || :verification in it.tags
        runs_loop && !tagged && push!(untagged_loops, Conventions.location(it))
        spawns_mpi && !(:mpi in it.tags) && push!(untagged_mpi, Conventions.location(it))
    end
    @test isempty(untagged_loops)
    @test isempty(untagged_mpi)
end

@testitem "conventions: per-commit simulations are the intended ones" tags=[:lint] setup=[Conventions] begin
    # The `:simulation` items that run on every commit. Adding a full simulation to CI means
    # editing this list, which is visible in review. Everything else with a time loop is
    # `:skipci` (verification) or `:mpi` (mpi job).
    allowed = Set([
        "NewtonKrylov throw maxiter",
        "Study: with a MultibodySetup",
        "contact!: two bodies of two points, one explicit step",
        "history-dependent model: the state is advanced once per step by a simulation",
        "material interface: a custom material runs a simulation",
        "process_each_export - serial and threads",
        "process_each_job: results, failed jobs and processing errors",
        "submit threads error handling",
        "submit!: all jobs succeed",
        "submit!: failing jobs are recorded and the others still run",
        "submit!: resuming skips the jobs a previous run completed",
        "submit!: the loud path prints and a failing job does not crash it",
        "symmetry: DynamicRelaxation",
        "symmetry: DynamicRelaxation, material variants",
        "symmetry: NewtonKrylov",
        "symmetry: NewtonKrylov, other materials",
        "symmetry: VelocityVerlet",
        "symmetry: VelocityVerlet, material variants",
        "uniform tension: DynamicRelaxation",
        "uniform tension: NewtonKrylov",
        "uniform tension: data boundary conditions",
    ])
    per_commit = Set(it.name for it in Conventions.all_items()
                     if :simulation in it.tags && !(:skipci in it.tags) && !(:mpi in it.tags))
    unexpected = setdiff(per_commit, allowed)
    missing_items = setdiff(allowed, per_commit)
    @test isempty(unexpected)
    @test isempty(missing_items)
end

@testitem "conventions: forbidden names" tags=[:lint] setup=[Conventions] begin
    # Files that still use forbidden names and are cleaned up one by one; the list shrinks to
    # nothing.
    pending = Set([
    ])
    violations = String[]
    for it in Conventions.all_items()
        isempty(it.lint_ok) || continue
        it.file in pending && continue
        for (name, reason) in Conventions.FORBIDDEN
            name in it.names || continue
            # random data is fine from the seeded generator `Fixtures.rng()`
            name in (:rand, :randn) && :rng in it.names && continue
            push!(violations, "$(Conventions.location(it)): $(name) — $(reason)")
        end
    end
    @test isempty(violations)
end

@testitem "conventions: unit test files mirror src" tags=[:lint] setup=[Conventions] begin
    # `test/unit/<dir>/test_<name>.jl` tests `src/<dir>/<name>.jl`: every directory below
    # `test/unit/` is a directory of `src/`, and every file in it names the src file it tests.
    # Files still to be renamed or moved are listed here and the list shrinks to nothing.
    pending = Set([])
    violations = String[]
    for dir in Conventions.unit_dirs()
        if !isdir(joinpath(Conventions.SRC_ROOT, dir))
            push!(violations, "unit/$(dir): no directory src/$(dir)")
            continue
        end
        for file in readdir(joinpath(Conventions.UNIT_ROOT, dir))
            endswith(file, ".jl") || continue
            rel = joinpath("unit", dir, file)
            rel in pending && continue
            m = match(r"^test_(.*)\.jl$", file)
            if m === nothing
                push!(violations, "$(rel): unit test files are named test_<srcfile>.jl")
                continue
            end
            src = joinpath(Conventions.SRC_ROOT, dir, m.captures[1] * ".jl")
            isfile(src) || push!(violations, "$(rel): no src file $(relpath(src, Conventions.SRC_ROOT))")
        end
    end
    @test isempty(violations)
end
