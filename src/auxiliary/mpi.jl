@inline mpi_comm() = MPI.COMM_WORLD
@inline mpi_rank() = MPI.Comm_rank(MPI.COMM_WORLD)
@inline mpi_nranks() = MPI.Comm_size(MPI.COMM_WORLD)
@inline mpi_run() = MPI_RUN[]
@inline mpi_chunk_id() = mpi_rank() + 1
@inline mpi_isroot() = MPI_ISROOT[]
@inline force_mpi_run!() = (MPI_RUN[] = true; MPI_RUN_FORCED[] = true; return nothing)

@inline function set_mpi_run!(b::Bool)
    MPI_RUN_FORCED[] && return nothing
    MPI_RUN[] = b
    return nothing
end

function init_mpi()
    if !MPI_INITIALIZED[]
        MPI.Init(finalize_atexit=true)
        MPI_INITIALIZED[] = true
    end
    MPI_ISROOT[] = mpi_rank() == 0
    set_mpi_run!(mpi_run_initial_check())
    return nothing
end

function mpi_run_initial_check()
    MPI_RUN_FORCED[] && return true
    haskey(ENV, "MPI_LOCALRANKID") && return true
    mpi_nranks() > 1 && return true
    nthreads() > 1 && return false
    for env_key in keys(ENV)
        if contains(env_key, "SLURM")
            msg = "will fall back to multithreading, due to "
            msg *= "`nranks` = `nthreads` = 1\n"
            msg *= "If you want to run MPI with only 1 rank, use `force_mpi_run!()`"
            @warn msg
            return false
        end
    end
    return false
end

function log_mpi_timers(options::ExportOptions)
    options.exportflag || return nothing
    file = joinpath(options.root, @sprintf("timers_rank_%d.log", mpi_rank()))
    open(file, "w+") do io
        show(IOContext(io, :displaysize => (24,150)), TO)
        write(io, "\n")
    end
    return nothing
end

function enable_mpi_timers()
    Core.eval(Peridynamics, :(TimerOutputs.enable_debug_timings(Peridynamics)))
    return nothing
end

function disable_mpi_timers()
    Core.eval(Peridynamics, :(TimerOutputs.disable_debug_timings(Peridynamics)))
    return nothing
end

# @inline function mpi_println(args...)
#     mpi_isroot() && println(args...)
#     return nothing
# end

# @inline function mpi_print(args...)
#     mpi_isroot() && print(args...)
#     return nothing
# end
