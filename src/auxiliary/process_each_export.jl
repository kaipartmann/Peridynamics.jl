const PROCESS_EACH_EXPORT_KWARGS = (:serial, :barrier)

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
    `serial=false` (parallel processing has automatic coordination). Use this when you need
    all ranks to synchronize after root-only file processing. (default: `false`)

# Returns
- `nothing` when `default_value === nothing` (legacy mode)
- `Vector{T}` where `T = typeof(default_value)` when result collection is enabled. Each
    element corresponds to the result from processing one export file. With MPI, all ranks
    receive the complete vector with results from all files.

# MPI Behavior
When running with MPI:
- `serial=true, barrier=false`: Processing runs only on root rank. Non-root ranks skip
    immediately. Use when synchronization is not needed or you'll add barriers manually.
- `serial=true, barrier=true`: Processing runs on root rank, then all ranks wait at a
    barrier. Use for the common case where all ranks need to wait for root's file I/O.
- `serial=false`: Processing is distributed across all MPI ranks with automatic
    coordination. Each rank processes a subset of files. The `barrier` parameter is ignored.
- With result collection enabled, `MPI.Allgatherv` is used to combine results from all ranks,
    ensuring every rank receives the complete results vector.

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
if mpi_isroot()
    println("Results from all files: ", results)
end
```

See also: [`process_each_job`](@ref), [`mpi_isroot`](@ref), [`mpi_barrier`](@ref)
"""
function process_each_export(f::F, vtk_path::AbstractString,
                             default_value=nothing; kwargs...) where {F<:Function}
    o = Dict{Symbol,Any}(kwargs)
    check_kwargs(o, PROCESS_EACH_EXPORT_KWARGS)
    check_process_function(f, default_value)
    serial, barrier = get_process_each_export_options(o)
    vtk_files = find_vtk_files(vtk_path)

    # Determine if we're collecting results
    collect_results = default_value !== nothing

    if serial
        results = @mpiroot process_each_export_serial(f, vtk_files, default_value)
        barrier && mpi_barrier()
        # In serial mode with result collection, broadcast results to all ranks
        if collect_results && mpi_run()
            results = broadcast_results(results, default_value, length(vtk_files))
        end
    elseif mpi_run()
        results = process_each_export_mpi(f, vtk_files, default_value)
    else
        results = process_each_export_threads(f, vtk_files, default_value)
    end

    return results
end

function process_each_export(f::F, job::Job, default_value=nothing;
                             kwargs...) where {F<:Function}
    return process_each_export(f, job.options.vtk, default_value; kwargs...)
end

function get_process_each_export_options(o::Dict{Symbol,Any})
    if haskey(o, :serial)
        serial::Bool = Bool(o[:serial])
    else
        serial = false
    end
    if haskey(o, :barrier)
        barrier::Bool = Bool(o[:barrier])
    else
        barrier = false
    end
    return serial, barrier
end

function check_process_function(f::F, default_value) where {F<:Function}
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

function process_step(f::F, ref_result::Dict{Symbol,T}, file::AbstractString,
                      file_id::Int, default_value) where {F<:Function,T}
    result = read_vtk(file)
    ret = try
        f(ref_result, result, file_id)
    catch err
        @error "something wrong while processing file $(basename(file))" error=err
        default_value
    end
    return ret
end

function process_each_export_serial(f::F, vtk_files::Vector{String},
                                    default_value) where {F<:Function}
    ref_result = read_vtk(first(vtk_files))
    p = Progress(length(vtk_files); dt=1, color=:normal, barlen=20,
                 enabled=progress_bars())

    collect_results = default_value !== nothing
    if collect_results
        T = typeof(default_value)
        results = Vector{T}(undef, length(vtk_files))
        for (file_id, file) in enumerate(vtk_files)
            results[file_id] = process_step(f, ref_result, file, file_id, default_value)
            next!(p)
        end
        finish!(p)
        return results
    else
        for (file_id, file) in enumerate(vtk_files)
            process_step(f, ref_result, file, file_id, nothing)
            next!(p)
        end
        finish!(p)
        return nothing
    end
end

function process_each_export_threads(f::F, vtk_files::Vector{String},
                                     default_value) where {F<:Function}
    ref_result = read_vtk(first(vtk_files))
    p = Progress(length(vtk_files); dt=1, color=:normal, barlen=20,
                 enabled=progress_bars())

    collect_results = default_value !== nothing
    if collect_results
        T = typeof(default_value)
        results = Vector{T}(undef, length(vtk_files))
        @threads :static for file_id in eachindex(vtk_files)
            results[file_id] = process_step(f, ref_result, vtk_files[file_id], file_id,
                                           default_value)
            next!(p)
        end
        finish!(p)
        return results
    else
        @threads :static for file_id in eachindex(vtk_files)
            process_step(f, ref_result, vtk_files[file_id], file_id, nothing)
            next!(p)
        end
        finish!(p)
        return nothing
    end
end

function process_each_export_mpi(f::F, vtk_files::Vector{String},
                                 default_value) where {F<:Function}
    ref_result = read_vtk(first(vtk_files))
    file_dist = distribute_equally(length(vtk_files), mpi_nranks())
    loc_file_ids = file_dist[mpi_chunk_id()]
    n_files = length(vtk_files)

    collect_results = default_value !== nothing

    if collect_results
        # Validate that result type is bits type for MPI
        T = typeof(default_value)
        if !isbitstype(T)
            msg = "result type must be a bits type for MPI compatibility!\n"
            msg *= "  Got type: $T\n"
            msg *= "  isbitstype($T) = $(isbitstype(T))\n"
            msg *= "Hint: Use NamedTuples of primitive types like Float64, Int, etc.\n"
            throw(ArgumentError(msg))
        end

        # Initialize results vector with default values
        results = fill(default_value, n_files)

        if mpi_isroot()
            p = Progress(length(vtk_files); dt=1, color=:normal, barlen=20,
                         enabled=progress_bars())
        end

        # Process local files
        for file_id in loc_file_ids
            results[file_id] = process_step(f, ref_result, vtk_files[file_id], file_id,
                                           default_value)
            mpi_isroot() && next!(p)
        end
        mpi_isroot() && finish!(p)

        # Gather results from all ranks
        results = gather_mpi_results(results, loc_file_ids, n_files)
        return results
    else
        if mpi_isroot()
            p = Progress(length(vtk_files); dt=1, color=:normal, barlen=20,
                         enabled=progress_bars())
        end
        for file_id in loc_file_ids
            process_step(f, ref_result, vtk_files[file_id], file_id, nothing)
            mpi_isroot() && next!(p)
        end
        mpi_isroot() && finish!(p)
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
    all_file_ids = Vector{Int}(undef, sum(counts))
    MPI.Allgatherv!(loc_file_ids, MPI.VBuffer(all_file_ids, counts), mpi_comm())

    # Reconstruct results in correct file order
    final_results = Vector{T}(undef, n_files)
    for (idx, file_id) in enumerate(all_file_ids)
        final_results[file_id] = recvbuf[idx]
    end

    return final_results
end

function broadcast_results(results, default_value, n_files::Int)
    # In serial mode, root has results, non-root ranks need to receive them
    if mpi_isroot()
        # Root broadcasts its results
        return results
    else
        # Non-root ranks receive the broadcast
        if default_value !== nothing
            ResultType = typeof(default_value)
            return Vector{ResultType}(undef, n_files)
        else
            return nothing
        end
    end
end
