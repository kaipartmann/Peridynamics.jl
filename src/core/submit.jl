
function submit(job::Job; kwargs...)
    o = Dict{Symbol,Any}(kwargs)
    check_kwargs(o, SUBMIT_KWARGS)
    quiet = get_submit_options(o)
    set_quiet(quiet)

    if mpi_sim()
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

function submit_mpi(::Job)
    error("functionality for MPI simulations not yet implemented!\n")
    return nothing
end

function submit_threads(job::Job, n_chunks::Int)
    point_decomp = PointDecomposition(job.spatial_setup, n_chunks)
    tdh = ThreadsDataHandler(job.spatial_setup, job.time_solver, point_decomp)
    init_time_solver!(job.time_solver, tdh)
    solve!(tdh, job.time_solver, job.options)
    return tdh
end
