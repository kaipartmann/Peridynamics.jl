using TestItemRunner

# The environment variable `PERIDYNAMICS_TESTS` selects which test items run, see
# `.github/workflows/CI.yml`:
#   "default" (unset): everything except the `:skipci`, `:mpi`, `:perf`, `:lint` and `:slow` items
#   "extras"         : only the `:mpi`, `:perf`, `:lint` and `:slow` items (without `:skipci`)
#   "all"            : everything except `:skipci`
# The `:verification` items are always `:skipci` and run deliberately with
#   julia -t 8 test/runtestitems.jl tag:verification
const SELECTION = get(ENV, "PERIDYNAMICS_TESTS", "default")
const EXTRA_TAGS = (:mpi, :perf, :lint, :slow)

function select(ti)
    :skipci in ti.tags && return false
    extra = any(in(ti.tags), EXTRA_TAGS)
    SELECTION == "default" && return !extra
    SELECTION == "extras" && return extra
    SELECTION == "all" && return true
    error("Unknown test selection PERIDYNAMICS_TESTS=$(SELECTION). " *
          "Use \"default\", \"extras\" or \"all\".")
end

@run_package_tests verbose=true filter=select
