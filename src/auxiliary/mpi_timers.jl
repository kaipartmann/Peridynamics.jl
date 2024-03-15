function log_mpi_timers(options::ExportOptions)
    options.exportflag || return nothing
    file = joinpath(options.root, @sprintf("timers_rank_%d.log", mpi_rank()))
    open(file, "w+") do io
        show(IOContext(io, :displaysize => (24,150)), TO)
        write(io, "\n")
    end
    return nothing
end

macro enable_mpi_timers()
    quote
        TimerOutputs.enable_debug_timings(Peridynamics)
        nothing
    end
end
