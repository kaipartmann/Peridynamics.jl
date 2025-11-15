const SUBMIT_KWARGS = (:quiet,)

"""
    submit(job::Job; quiet=false)

Run the simulation by submitting the job.

# Arguments
- `job::Job`: Job that contains all defined parameters and conditions.

# Keywords
- `quiet::Bool`: If `true`, no outputs are printed in the terminal. (default: `false`)
"""
function submit(job::Job; kwargs...)
    o = Dict{Symbol,Any}(kwargs)
    check_kwargs(o, SUBMIT_KWARGS)
    quiet = get_submit_options(o)
    set_quiet!(quiet)

    if mpi_run()
        dh = submit_mpi(job)
    else
        dh = submit_threads(job, nthreads())
    end

    return dh
end

function get_submit_options(o::Dict{Symbol,Any})
    local quiet::Bool
    if haskey(o, :quiet)
        quiet = Bool(o[:quiet])
    else
        quiet = false
    end
    return quiet
end

function submit_mpi(job::Job)
    dh = try
        _submit_mpi(job)
    catch err
        timestamp = Dates.format(Dates.now(), "yyyy-mm-dd, HH:MM:SS")
        msg = "\nERROR: submit failed with error at $(timestamp) in rank $(mpi_rank())!\n"
        msg *= sprint(showerror, err, catch_backtrace())
        msg *= "\n"
        add_to_logfile(job.options, msg)
        # only rethrowing NaNErrors, which are synchronized across all ranks
        isa(err, NaNError) && rethrow(err)
        # for other errors, abort MPI to avoid deadlocks (if more than one rank exists)
        mpi_nranks() > 1 && MPI.Abort(mpi_comm(), 1)
        # if it is happening only on a single rank, we can rethrow safely
        rethrow(err)
    end
    return dh
end

function _submit_mpi(job::Job)
    timeit_debug_enabled() && reset_timer!(TO)
    simulation_duration = @elapsed begin
        @timeit_debug TO "initialization" begin
            init_logs(job.options)
            log_spatial_setup(job.options, job.spatial_setup)
            log_create_data_handler_start()
            dh = mpi_data_handler(job.spatial_setup, job.time_solver)
            init_time_solver!(job.time_solver, dh)
            initialize!(dh, job.time_solver)
            log_create_data_handler_end()
            log_data_handler(job.options, dh)
            log_job_options(job.options)
            log_timesolver(job.options, job.time_solver)
        end
        @timeit_debug TO "solve!" begin
            solve!(dh, job.time_solver, job.options)
        end
    end
    log_simulation_duration(job.options, simulation_duration)
    if timeit_debug_enabled()
        TimerOutputs.complement!(TO)
        log_mpi_timers(job.options)
    end
    return dh
end

function submit_threads(job::Job, n_threads::Int)
    dh = try
        _submit_threads(job, n_threads)
    catch err
        timestamp = Dates.format(Dates.now(), "yyyy-mm-dd, HH:MM:SS")
        msg = "\nERROR: submit failed with error at $(timestamp)!\n"
        msg *= sprint(showerror, err, catch_backtrace())
        msg *= "\n"
        add_to_logfile(job.options, msg)
        rethrow(err) # rethrowing all errors, as no deadlocks can occur here
    end
    return dh
end

function _submit_threads(job::Job, n_chunks::Int)
    simulation_duration = @elapsed begin
        init_logs(job.options)
        log_spatial_setup(job.options, job.spatial_setup)
        log_create_data_handler_start()
        dh = threads_data_handler(job.spatial_setup, job.time_solver, n_chunks)
        init_time_solver!(job.time_solver, dh)
        initialize!(dh, job.time_solver)
        log_create_data_handler_end()
        log_data_handler(job.options, dh)
        log_job_options(job.options)
        log_timesolver(job.options, job.time_solver)
        solve!(dh, job.time_solver, job.options)
    end
    log_simulation_duration(job.options, simulation_duration)
    return dh
end
