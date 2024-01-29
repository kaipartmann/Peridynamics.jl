
function submit(job::J; nchunks::Int=nthreads()) where {J<:AbstractJob}
    if mpi_sim()
        ret = submit_mpi(job)
    else
        ret = submit_threads(job, nchunks)
    end
    return ret
end

"""
MAIN FUNCTION THREADS
"""
function submit_threads(job::J, n_chunks::Int) where {J<:AbstractJob}
    tdh = init_threads_data_handler(job, n_chunks)
end
