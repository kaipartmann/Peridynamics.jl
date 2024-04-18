const SUBMIT_KWARGS = (:quiet,)

"""
    submit(job::Job; kwargs...)

submits the job to start calculations

# Arguments

- `job::Job`: job that contains all defined parameters and conditions

# Keywords

- `:quiet::Bool`: if `:quiet=true`, no outputs are printed in the terminal

# Throws

- error if keyword is not allowed

# Example

```julia-repl
julia> job = Job(b, vv;
           path=joinpath(@__DIR__, "results", "mode_I"),
           fields=(:displacement, :damage))
julia> submit(job)
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
    @timeit_debug TO "initialization" begin
        init_logs(job.options)
        point_decomp = PointDecomposition(job.spatial_setup, mpi_nranks())
        mdh = MPIDataHandler(job.spatial_setup, job.time_solver, point_decomp)
        init_time_solver!(job.time_solver, mdh)
        init_data_handler!(mdh, job.time_solver)
    end
    @timeit_debug TO "solve!" begin
        solve!(mdh, job.time_solver, job.options)
    end
    if timeit_debug_enabled()
        TimerOutputs.complement!(TO)
        log_mpi_timers(job.options)
    end
    return mdh
end

function submit_threads(job::Job, n_chunks::Int)
    init_logs(job.options)
    point_decomp = PointDecomposition(job.spatial_setup, n_chunks)
    tdh = ThreadsDataHandler(job.spatial_setup, job.time_solver, point_decomp)
    init_time_solver!(job.time_solver, tdh)
    init_data_handler!(tdh, job.time_solver)
    solve!(tdh, job.time_solver, job.options)
    return tdh
end
