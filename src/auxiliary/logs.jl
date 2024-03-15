is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")

@inline progress_enabled() = !is_logging(stderr) && !quiet()

function init_logs(options::ExportOptions)
    options.exportflag || return nothing
    mkpath(options.vtk)
    init_logfile(options)
    return nothing
end

function init_logfile(options::ExportOptions)
    open(options.logfile, "w+") do io
        write(io, "--- LOGFILE ---\n")
    end
    return nothing
end
