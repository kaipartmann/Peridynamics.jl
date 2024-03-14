is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")

@inline progress_enabled() = !is_logging(stderr) && !quiet()

function log_timeroutputs(options::ExportOptions)
    options.exportflag || return nothing
    open(options.logfile, "a+") do io
        write(io, "\n")
        show(IOContext(io, :displaysize => (24,150)), TO)
    end
    return nothing
end
