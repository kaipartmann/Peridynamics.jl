using TestItemRunner

# The environment variable `PERIDYNAMICS_TESTS` selects which test items run, see
# `.github/workflows/CI.yml`:
#   "default" (unset): everything except the `:skipci`, `:mpi`, `:perf`, `:lint` and `:slow` items
#   "mpi"            : only the `:mpi` items (without `:skipci`)
#   "extras"         : only the `:perf`, `:lint` and `:slow` items (without `:skipci`)
#   "all"            : everything except `:skipci`
# The `:verification` items are always `:skipci` and run deliberately with
#   julia -t 8 test/runtestitems.jl tag:verification
const SELECTION = get(ENV, "PERIDYNAMICS_TESTS", "default")
const EXTRA_TAGS = (:perf, :lint, :slow)

function select(ti)
    :skipci in ti.tags && return false
    mpi = :mpi in ti.tags
    extra = any(in(ti.tags), EXTRA_TAGS)
    SELECTION == "default" && return !(mpi || extra)
    SELECTION == "mpi" && return mpi
    SELECTION == "extras" && return extra && !mpi
    SELECTION == "all" && return true
    error("Unknown test selection PERIDYNAMICS_TESTS=$(SELECTION). " *
          "Use \"default\", \"mpi\", \"extras\" or \"all\".")
end

@run_package_tests verbose=true filter=select
