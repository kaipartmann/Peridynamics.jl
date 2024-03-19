function mpi_run_check()
    haskey(ENV, "MPI_LOCALRANKID") && return true
    mpi_nranks() > 1 && return true
    if mpi_nranks() == 1 && nthreads() == 1
        msg = "unsure if this is a MPI or multithreading run, due to "
        msg *= "`nranks` == `nthreads` == 1\n"
        msg *= "The fallback assumption for this case is: MPI"
        @warn msg
        return true
    end
    return false
end
