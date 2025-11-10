const PROCESS_EACH_EXPORT_KWARGS = (:serial, :barrier)

"""
    process_each_export(f, vtk_path; kwargs...)
    process_each_export(f, job; kwargs...)

A function for postprocessing every exported file. This function works with multithreading
and MPI and determines the backend exactly like the [`submit`](@ref) function.

# Arguments
- `f::Function`: The processing function with signature `f(r0, r, id)`.
    - `r0`: The results of [`read_vtk`](@ref) for the exported file of the reference
        results.
    - `r`: The results of [`read_vtk`](@ref) for a time step.
    - `id::Ind`: An ID indicating the number of the exported file (counted from 1, starting
        with the reference file).
- `vtk_path::AbstractString`: A path that should contain the export results of a simulation.
- `job::Job`: A job object. The path of the VTK files will then be processed from the
    job options.

# Keywords
- `serial::Bool`: If `true`, all results will be processed in the correct order of the time
    steps and on a single thread, cf. the MPI root rank. (default: `false`)
- `barrier::Bool`: If `true` and `serial=true`, adds an MPI barrier after processing
    completes, ensuring all ranks wait for the root rank to finish. This is ignored when
    `serial=false` (parallel processing has automatic coordination). Use this when you need
    all ranks to synchronize after root-only file processing. (default: `false`)

# MPI Behavior
When running with MPI:
- `serial=true, barrier=false`: Processing runs only on root rank. Non-root ranks skip
    immediately. Use when synchronization is not needed or you'll add barriers manually.
- `serial=true, barrier=true`: Processing runs on root rank, then all ranks wait at a
    barrier. Use for the common case where all ranks need to wait for root's file I/O.
- `serial=false`: Processing is distributed across all MPI ranks with automatic
    coordination. Each rank processes a subset of files. The `barrier` parameter is ignored.

# Examples
```julia
# Root processes files, all ranks wait
process_each_export(job; serial=true, barrier=true) do r0, r, id
    # File operations on root
end
# All ranks synchronized here

# Root processes files, no waiting
process_each_export(job; serial=true, barrier=false) do r0, r, id
    # File operations on root
end
# Non-root ranks continue immediately

# All ranks process files in parallel
process_each_export(job) do r0, r, id
    # File operations on all ranks
end
```

See also: [`process_each_job`](@ref), [`mpi_isroot`](@ref), [`mpi_barrier`](@ref)
"""
function process_each_export(f::F, vtk_path::AbstractString; kwargs...) where {F<:Function}
    o = Dict{Symbol,Any}(kwargs)
    check_kwargs(o, PROCESS_EACH_EXPORT_KWARGS)
    check_process_function(f)
    serial, barrier = get_process_each_export_options(o)
    vtk_files = find_vtk_files(vtk_path)
    if serial
        @mpiroot process_each_export_serial(f, vtk_files)
        barrier && mpi_barrier()
    elseif mpi_run()
        process_each_export_mpi(f, vtk_files)
    else
        process_each_export_threads(f, vtk_files)
    end
    return nothing
end

function process_each_export(f::F, job::Job; kwargs...) where {F<:Function}
    process_each_export(f, job.options.vtk; kwargs...)
    return nothing
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

function process_step(f::F, ref_result::Dict{Symbol,T}, file::AbstractString,
                      file_id::Int) where {F<:Function,T}
    result = read_vtk(file)
    try
        f(ref_result, result, file_id)
    catch err
        @error "something wrong while processing file $(basename(file))" error=err
    end
    return nothing
end

function process_each_export_serial(f::F, vtk_files::Vector{String}) where {F<:Function}
    ref_result = read_vtk(first(vtk_files))
    p = Progress(length(vtk_files); dt=1, color=:normal, barlen=20,
                 enabled=progress_bars())
    for (file_id, file) in enumerate(vtk_files)
        process_step(f, ref_result, file, file_id)
        next!(p)
    end
    finish!(p)
    return nothing
end

function process_each_export_threads(f::F, vtk_files::Vector{String}) where {F<:Function}
    ref_result = read_vtk(first(vtk_files))
    p = Progress(length(vtk_files); dt=1, color=:normal, barlen=20,
                 enabled=progress_bars())
    @threads :static for file_id in eachindex(vtk_files)
        process_step(f, ref_result, vtk_files[file_id], file_id)
        next!(p)
    end
    finish!(p)
    return nothing
end

function process_each_export_mpi(f::F, vtk_files::Vector{String}) where {F<:Function}
    ref_result = read_vtk(first(vtk_files))
    file_dist = distribute_equally(length(vtk_files), mpi_nranks())
    loc_file_ids = file_dist[mpi_chunk_id()]
    if mpi_isroot()
        p = Progress(length(vtk_files); dt=1, color=:normal, barlen=20,
                     enabled=progress_bars())
    end
    for file_id in loc_file_ids
        process_step(f, ref_result, vtk_files[file_id], file_id)
        mpi_isroot() && next!(p)
    end
    mpi_isroot() && finish!(p)
    return nothing
end
