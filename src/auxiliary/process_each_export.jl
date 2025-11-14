"""
    process_each_export(f, vtk_path, default_value=nothing; kwargs...)
    process_each_export(f, job, default_value=nothing; kwargs...)

A function for postprocessing every exported file. This function works with multithreading
and MPI and determines the backend exactly like the [`submit`](@ref) function.

# Arguments
- `f::Function`: The processing function with signature `f(r0, r, id)` that returns a result
    to be collected (when `default_value !== nothing`) or returns `nothing` (legacy mode).
    - `r0`: The results of [`read_vtk`](@ref) for the exported file of the reference
        results.
    - `r`: The results of [`read_vtk`](@ref) for a time step.
    - `id::Int`: An ID indicating the number of the exported file (counted from 1, starting
        with the reference file).
- `vtk_path::AbstractString`: A path that should contain the export results of a simulation.
- `job::Job`: A job object. The path of the VTK files will then be processed from the
    job options.
- `default_value`: Optional default value for result collection. When provided (not `nothing`),
    the function returns a vector of results from processing each export file. The type of
    `default_value` determines the element type of the returned vector. For MPI compatibility,
    the returned result type must satisfy `isbitstype(result) == true`.
    When `default_value === nothing` (default), the function returns `nothing` and behaves
    like the original version without result collection.

# Keywords
- `serial::Bool`: If `true`, all results will be processed in the correct order of the time
    steps and on a single thread, cf. the MPI root rank. (default: `false`)
- `barrier::Bool`: If `true` and `serial=true`, adds an MPI barrier after processing
    completes, ensuring all ranks wait for the root rank to finish. This is ignored when
    `serial=false` (parallel processing has automatic coordination) or when `only_root=true`.
    Use this when you need all ranks to synchronize after root-only file processing.
    (default: `false`)
- `only_root::Bool`: If `true`, processing runs only on the root rank and results are NOT
    broadcast to other ranks. This prevents MPI deadlocks when calling `process_each_export`
    inside `process_each_job` with `only_root=true`. Non-root ranks return `nothing`.
    When `false`, results are broadcast to all ranks (if result collection is enabled).
    (default: `false`)

# Returns
- `nothing` when `default_value === nothing` (legacy mode)
- `Vector{T}` where `T = typeof(default_value)` when result collection is enabled. Each
    element corresponds to the result from processing one export file. With MPI, all ranks
    receive the complete vector with results from all files.

# MPI Behavior
When running with MPI:
- `serial=true, barrier=false, only_root=false`: Processing runs only on root rank. Results
    are broadcast to all ranks. Non-root ranks wait for broadcast but don't process files.
- `serial=true, barrier=true, only_root=false`: Processing runs on root rank, results are
    broadcast to all ranks, then all ranks wait at a barrier.
- `only_root=true`: Processing runs only on root rank. Results are NOT broadcast - non-root
    ranks return `nothing` immediately. Use this to avoid MPI deadlocks when calling
    `process_each_export` inside `process_each_job(...; only_root=true)`. The `barrier`
    parameter is ignored.
- `serial=false`: Processing is distributed across all MPI ranks with automatic coordination.
    Each rank processes a subset of files. The `barrier` and `only_root` parameters are
    ignored.
- With result collection enabled and `only_root=false`, results are gathered/broadcast so
    all ranks receive the complete results vector. With `only_root=true`, only the root rank
    gets results.

!!! warning "MPI result type requirements"
    When using result collection with MPI (`default_value !== nothing`), the returned result
    type must be a bits type (`isbitstype(result) == true`). This typically means primitive
    types or `NamedTuple`s of primitive types. Non-bits types like strings or arrays will
    cause an error.

# Examples
```julia
# Legacy mode - no result collection
process_each_export(job; serial=true, barrier=true) do r0, r, id
    # File operations on root
    open("result_\$id.txt", "w") do io
        write(io, string(maximum(r[:displacement])))
    end
end

# Result collection mode
default_value = (; max_disp=NaN, avg_ux=NaN)
results = process_each_export(job, default_value) do r0, r, id
    max_disp = maximum(r[:displacement])
    ux = @view r[:displacement][1, :]
    avg_ux = sum(ux) / length(ux)
    return (; max_disp, avg_ux)
end
# results is a Vector{NamedTuple{(:max_disp, :avg_ux), Tuple{Float64, Float64}}}

# With MPI, all ranks get the complete results
@mpiroot println("Results from all files: ", results)

# Using inside process_each_job to avoid deadlocks
process_each_job(jobs; only_root=true) do job, job_id
    submit(job)
    # Use only_root=true here to prevent MPI deadlock
    results = process_each_export(job, default_value; only_root=true) do r0, r, id
        # Process export files...
        return (; max_disp=maximum(r[:displacement]))
    end
    # Only root rank has results here, non-root ranks get nothing
    return results
end
```

See also: [`process_each_job`](@ref), [`mpi_isroot`](@ref), [`mpi_barrier`](@ref)
"""
function process_each_export(f::F, vtk_path::AbstractString, default_ret=nothing;
                             serial::Bool=false, barrier::Bool=false,
                             only_root::Bool=false) where {F<:Function}
    check_process_function(f)
    vtk_files = find_vtk_files(vtk_path)

    if serial || only_root
        results = @mpiroot process_each_export_serial(f, vtk_files, default_ret)
        barrier && !only_root && mpi_barrier()
        # In serial mode with result collection, broadcast results to all ranks
        if default_ret !== nothing && mpi_run() && !only_root
            results = broadcast_results(results, default_ret, length(vtk_files))
        end
    elseif mpi_run()
        results = process_each_export_mpi(f, vtk_files, default_ret)
    else
        results = process_each_export_threads(f, vtk_files, default_ret)
    end

    return results
end

function process_each_export(f::F, job::Job, default_ret=nothing;
                             kwargs...) where {F<:Function}
    return process_each_export(f, job.options.vtk, default_ret; kwargs...)
end

function check_process_function(f::F) where {F<:Function}
    func_method = get_method_of_function(f)
    args = get_argument_names_of_function(func_method)
    if length(args) != 3
        msg = "wrong arguments for processing function!\n"
        msg *= "The processing function needs 3 args:\n"
        msg *= "  - ref_result::Dict{Symbol,VecOrMat}: results of the initial export\n"
        msg *= "  - result::Dict{Symbol,VecOrMat}: results of the processed time step\n"
        msg *= "  - file_id::Int: number / id of the processed export\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function find_vtk_files(path::AbstractString)
    isdir(path) || throw(ArgumentError("invalid path $path specified!\n"))
    all_files = readdir(path; join=true)
    pvtu_files = filter(x -> endswith(x, ".pvtu"), all_files)
    isempty(pvtu_files) && throw(ArgumentError("no pvtu-files in path: $path\n"))
    sort_by_time_step!(pvtu_files)
    return pvtu_files
end

function extract_time_step_from_filename(filename::AbstractString)
    filename_base = first(splitext(basename(filename)))
    time_step_str = last(split(filename_base, "_"))
    time_step = parse(Int, time_step_str)
    #TODO
    # if startswith(filename_base, "timestep")
    #     body_name = ""
    # else
    #     body_name = first(split(filename_base, "_timestep"))
    # end
    # return body_name, time_step
    return time_step
end

function sort_by_time_step!(vtk_files::Vector{<:AbstractString})
    idxs = sortperm(extract_time_step_from_filename.(vtk_files))
    vtk_files .= vtk_files[idxs]
    return nothing
end

function process_step(f::F, ref_result::Dict{Symbol,T}, vtk_file::AbstractString, id::Int,
                      default_ret) where {F<:Function,T}
    error_occurred = false
    ret = try
        f(ref_result, read_vtk(vtk_file), id)
    catch err
        error_occurred = true
        write_processing_error_log(vtk_file, id, err, catch_backtrace())
        default_ret
    end
    return ret, error_occurred
end

function write_processing_error_log(vtk_file::AbstractString, id::Int, err, bt)
    err_dir = joinpath(dirname(vtk_file), "vtk_process_errors")
    isdir(err_dir) || mkpath(err_dir)
    filename = @sprintf("proc_error_file_%d.log", id)
    logfile = joinpath(err_dir, filename)
    err_msg = "ERROR: vtk processing failed with error!\n"
    err_msg *= "File: $(basename(vtk_file))\n"
    err_msg *= "File ID: $(id)\n\n"
    err_msg *= sprint(showerror, err, bt)
    open(logfile, "w+") do io
        write(io, get_logfile_head())
        write(io, peridynamics_banner(color=false))
        write(io, err_msg)
    end
    return nothing
end

function check_proc_errors(n_errors::Int)
    mpi_isroot() || return nothing
    if n_errors > 0
        err_dir = joinpath("vtk", "vtk_process_errors")
        err_msg = @sprintf("%d error(s) occurred during export processing!\n", n_errors)
        err_msg *= "       See detailed error logs in: $(err_dir)\n"
        printstyled(stderr, "\nERROR: "; color=:red, bold=true)
        print(stderr, err_msg)
    end
    return nothing
end
check_proc_errors(::Nothing) = nothing

function process_each_export_serial(f::F, vtk_files::Vector{String},
                                    default_ret::T) where {F<:Function,T}
    ref_result = read_vtk(first(vtk_files))
    p = Progress(length(vtk_files); dt=1, color=:normal, barlen=20,
                 enabled=progress_bars())
    if default_ret !== nothing
        error_count = 0
        rets = Vector{T}(undef, length(vtk_files))
        for (id, vtk_file) in enumerate(vtk_files)
            ret, error_occurred = process_step(f, ref_result, vtk_file, id, default_ret)
            rets[id] = ret
            error_occurred && (error_count += 1)
            next!(p)
        end
        finish!(p)
        check_proc_errors(error_count)
        return rets
    else
        error_count = 0
        for (id, vtk_file) in enumerate(vtk_files)
            _, error_occurred = process_step(f, ref_result, vtk_file, id, nothing)
            error_occurred && (error_count += 1)
            next!(p)
        end
        finish!(p)
        check_proc_errors(error_count)
        return nothing
    end
end

function process_each_export_threads(f::F, vtk_files::Vector{String},
                                     default_ret::T) where {F<:Function,T}
    ref_result = read_vtk(first(vtk_files))
    p = Progress(length(vtk_files); dt=1, color=:normal, barlen=20,
                 enabled=progress_bars())
    if default_ret !== nothing
        error_count = Atomic{Int}(0) # thread-safe error counting using atomic
        rets = Vector{T}(undef, length(vtk_files))
        @threads for id in eachindex(vtk_files)
            vtk_file = vtk_files[id]
            ret, error_occurred = process_step(f, ref_result, vtk_file, id, default_ret)
            rets[id] = ret
            error_occurred && atomic_add!(error_count, 1)
            next!(p)
        end
        finish!(p)
        check_proc_errors(error_count[])
        return rets
    else
        error_count = Atomic{Int}(0) # thread-safe error counting using atomic
        @threads for id in eachindex(vtk_files)
            vtk_file = vtk_files[id]
            _, error_occurred = process_step(f, ref_result, vtk_file, id, nothing)
            error_occurred && atomic_add!(error_count, 1)
            next!(p)
        end
        finish!(p)
        check_proc_errors(error_count[])
        return nothing
    end
end

function process_each_export_mpi(f::F, vtk_files::Vector{String},
                                 default_ret) where {F<:Function}
    ref_result = read_vtk(first(vtk_files))
    file_dist = distribute_equally(length(vtk_files), mpi_nranks())
    loc_file_ids = file_dist[mpi_chunk_id()]
    n_files = length(vtk_files)
    if default_ret !== nothing
        loc_error_count = 0 # count errors locally on each rank
        # Validate that result type is bits type for MPI
        check_default_value_isbitstype(default_ret)
        # Initialize results vector with default values
        rets = fill(default_ret, n_files)
        if mpi_isroot()
            p = Progress(length(loc_file_ids); dt=1, color=:normal, barlen=20,
                         enabled=progress_bars())
        end
        # Process local files
        for id in loc_file_ids
            vtk_file = vtk_files[id]
            ret, error_occurred = process_step(f, ref_result, vtk_file, id, default_ret)
            rets[id] = ret
            error_occurred && (loc_error_count += 1)
            mpi_isroot() && next!(p)
        end
        mpi_isroot() && finish!(p)
        # Gather results from all ranks
        rets = gather_mpi_results(rets, loc_file_ids, n_files)
        # Gather error counts and print summary only on root
        total_errors = MPI.Reduce(loc_error_count, +, mpi_comm(); root=0)
        check_proc_errors(total_errors)
        mpi_barrier()
        return rets
    else
        loc_error_count = 0 # count errors locally on each rank
        if mpi_isroot()
            p = Progress(length(loc_file_ids); dt=1, color=:normal, barlen=20,
                         enabled=progress_bars())
        end
        for id in loc_file_ids
            vtk_file = vtk_files[id]
            _, error_occurred = process_step(f, ref_result, vtk_file, id, nothing)
            error_occurred && (loc_error_count += 1)
            mpi_isroot() && next!(p)
        end
        mpi_isroot() && finish!(p)
        # Gather error counts and print summary only on root
        total_errors = MPI.Reduce(loc_error_count, +, mpi_comm(); root=0)
        check_proc_errors(total_errors)
        mpi_barrier()
        return nothing
    end
end

function gather_mpi_results(results::Vector{T}, loc_file_ids::AbstractVector{Int},
                            n_files::Int) where {T}
    # Create send buffer with only local results
    n_loc = length(loc_file_ids)
    sendbuf = results[loc_file_ids]

    # Gather counts from all ranks
    counts = MPI.Allgather(n_loc, mpi_comm())

    # Create receive buffer
    recvbuf = Vector{T}(undef, sum(counts))

    # Perform Allgatherv
    MPI.Allgatherv!(sendbuf, MPI.VBuffer(recvbuf, counts), mpi_comm())

    # Gather file IDs from all ranks to reconstruct correct order
    # Convert to Vector{Int} in case loc_file_ids is a range
    loc_file_ids_vec = collect(Int, loc_file_ids)
    all_file_ids = Vector{Int}(undef, sum(counts))
    MPI.Allgatherv!(loc_file_ids_vec, MPI.VBuffer(all_file_ids, counts), mpi_comm())

    # Reconstruct results in correct file order
    final_results = Vector{T}(undef, n_files)
    for (idx, file_id) in enumerate(all_file_ids)
        final_results[file_id] = recvbuf[idx]
    end

    return final_results
end

function broadcast_results(results, default_value, n_files::Int)
    # Type must be bits type for MPI broadcast
    check_default_value_isbitstype(default_value)

    # Create or use existing buffer
    if mpi_isroot()
        # Root already has results from processing
        buffer = results
    else
        # Non-root ranks create buffer to receive results
        buffer = Vector{typeof(default_value)}(undef, n_files)
    end

    # Broadcast from root (rank 0) to all ranks
    MPI.Bcast!(buffer, mpi_comm(); root=0)

    return buffer
end

function check_default_value_isbitstype(::T) where {T}
    if !isbitstype(T)
        msg = "result type must be a bits type for MPI compatibility!\n"
        msg *= "  Got type: $T\n"
        msg *= "  isbitstype($T) = $(isbitstype(T))\n"
        msg *= "Hint: Use NamedTuples of primitive types like Float64, Int, etc.\n"
        throw(ArgumentError(msg))
    end
    return nothing
end
