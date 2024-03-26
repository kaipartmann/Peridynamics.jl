@inline quiet() = QUIET[]
@inline set_quiet!(b::Bool) = (QUIET[] = b; return nothing)

is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")

@inline progress_enabled() = !is_logging(stderr) && !quiet()

function init_logs(options::ExportOptions)
    options.exportflag || return nothing
    mpi_isroot() || return nothing
    mkpath(options.vtk)
    init_logfile(options)
    return nothing
end

function init_logfile(options::ExportOptions)
    open(options.logfile, "w+") do io
        write(io, "--- LOGFILE ---\n")
        if mpi_run()
            write(io, @sprintf("MPI SIMULATION WITH %d RANKS\n", mpi_nranks()))
        else
            write(io, @sprintf("MULTITHREADING SIMULATION WITH %d THREADS\n", nthreads()))
        end
    end
    return nothing
end
