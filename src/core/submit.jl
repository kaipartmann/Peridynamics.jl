const SUBMIT_KWARGS = (:quiet,)

"""
    submit(job::Job; kwargs...)

Submits the job to start calculations

# Arguments

- `job::Job`: Job that contains all defined parameters and conditions

# Keywords

- `quiet::Bool`: If `:quiet=true`, no outputs are printed in the terminal.

# Throws

- Error if keyword is not allowed

# Example

```julia-repl
julia> job = Job(b, vv;
           path=joinpath(@__DIR__, "results", "mode_I"),
           fields=(:displacement, :damage));

julia> submit(job);
solve... 100%|████████████████████| Time: 0:00:11
```
"""
function submit(job::Job; kwargs...)
    o = Dict{Symbol,Any}(kwargs)
    check_kwargs(o, SUBMIT_KWARGS)
    quiet = get_submit_options(o)
    set_quiet!(quiet)

    if mpi_run()
        ret = submit_mpi(job)
    else
        ret = submit_threads(job, nthreads())
    end

    return ret
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

function submit_threads(job::Job, n_chunks::Int)
    simulation_duration = @elapsed begin
        init_logs(job.options)
        log_spatial_setup(job.options, job.spatial_setup)
        log_create_data_handler_start()
        dh = threads_data_handler(job.spatial_setup, job.time_solver, n_chunks)
        init_time_solver!(job.time_solver, dh)
        initialize!(dh, job.time_solver)
        log_create_data_handler_end()
        log_data_handler(job.options, dh)
        log_timesolver(job.options, job.time_solver)
        solve!(dh, job.time_solver, job.options)
    end
    log_simulation_duration(job.options, simulation_duration)
    return dh
end
