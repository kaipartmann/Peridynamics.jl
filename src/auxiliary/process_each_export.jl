function process_each_export(f::F, vtk_path::AbstractString; kwargs...) where {F<:Function}
    o = Dict{Symbol,Any}(kwargs)
    check_kwargs(o, PROCESS_EACH_EXPORT_KWARGS)
    serial = get_process_each_export_options(o)
    vtk_files = find_vtk_files(vtk_path)
    if serial
        process_each_export_serial(f, vtk_files)
    elseif mpi_sim()
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
    return serial
end

function find_vtk_files(path::AbstractString)
    isdir(path) || throw(ArgumentError("invalid path $path specified!\n"))
    all_files = readdir(path; join=true)
    pvtu_files = filter(x -> endswith(x, ".pvtu"), all_files)
    isempty(pvtu_files) && throw(ArgumentError("no pvtu-files in path: $path\n"))
    return pvtu_files
end

function process_step(f::F, file::AbstractString, file_id::Int) where {F<:Function}
    result = read_vtk(file)
    try
        f(result, file_id)
    catch err
        @error "something wrong while processing file $(basename(file))" error=err
    end
    return nothing
end

function process_each_export_serial(f::F, vtk_files::Vector{String}) where {F<:Function}
    p = Progress(length(vtk_files); dt=1, color=:normal, barlen=20,
                 enabled=progress_enabled())
    for (file_id, file) in enumerate(vtk_files)
        process_step(f, file, file_id)
        next!(p)
    end
    finish!(p)
    return nothing
end

function process_each_export_threads(f::F, vtk_files::Vector{String}) where {F<:Function}
    p = Progress(length(vtk_files); dt=1, color=:normal, barlen=20,
                 enabled=progress_enabled())
    @threads for file_id in eachindex(vtk_files)
        process_step(f, vtk_files[file_id], file_id)
        next!(p)
    end
    finish!(p)
    return nothing
end

function process_each_export_mpi(f::F, vtk_files::Vector{String}) where {F<:Function}
    file_dist = distribute_equally(length(vtk_files), mpi_nranks())
    loc_file_ids = file_dist[mpi_rank()+1]
    if mpi_rank() == 0
        p = Progress(length(vtk_files); dt=1, color=:normal, barlen=20,
                     enabled=progress_enabled())
    end
    for file_id in loc_file_ids
        process_step(f, vtk_files[file_id], file_id)
        mpi_rank() == 0 && next!(p)
    end
    mpi_rank() == 0 && finish!(p)
    return nothing
end
