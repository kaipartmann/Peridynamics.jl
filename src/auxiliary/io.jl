const ExportField = Tuple{Symbol,DataType}

struct ExportOptions{N}
    exportflag::Bool
    root::String
    vtk::String
    logfile::String
    freq::Int
    fields::NTuple{N,ExportField}

    function ExportOptions(root::String, freq::Int, fields::NTuple{N,ExportField}) where {N}
        if isempty(root)
            return new{0}(false, "", "", "", 0, NTuple{0,ExportField}())
        end
        vtk = joinpath(root, "vtk")
        logfile = joinpath(root, "logfile.log")
        return new{N}(true, root, vtk, logfile, freq, fields)
    end
end

function get_export_options(::Type{S}, o::Dict{Symbol,Any}) where {S<:AbstractStorage}
    local root::String
    local freq::Int

    if haskey(o, :path) && haskey(o, :freq)
        root = string(o[:path])
        freq = Int(o[:freq])
    elseif haskey(o, :path) && !haskey(o, :freq)
        root = string(o[:path])
        freq = 10
    elseif !haskey(o, :path) && haskey(o, :freq)
        msg = "if `freq` is spedified, the keyword `path` is also needed!\n"
        throw(ArgumentError(msg))
    else
        root = ""
        freq = 0
    end
    freq < 0 && throw(ArgumentError("`freq` should be larger than zero!\n"))

    fields = get_export_fields(S, o)

    return ExportOptions(root, freq, fields)
end

function get_export_fields(::Type{S}, o::Dict{Symbol,Any}) where {S}
    export_fieldnames = get_export_fieldnames(o)
    storage_fieldnames = fieldnames(S)
    storage_fieldtypes = fieldtypes(S)

    _export_fields = Vector{ExportField}()

    for name in export_fieldnames
        idx = findfirst(x -> x === name, storage_fieldnames)
        isnothing(idx) && unknown_fieldname_error(S, name)
        type = storage_fieldtypes[idx]
        push!(_export_fields, (name, type))
    end

    export_fields = Tuple(_export_fields)

    return export_fields
end

function unknown_fieldname_error(::Type{S}, name::Symbol) where {S<:AbstractStorage}
    msg = "unknown field $(name) specified for export!\n"
    msg *= "Allowed fields for $S:\n"
    for allowed_name in fieldnames(S)
        msg *= "  - $allowed_name\n"
    end
    throw(ArgumentError(msg))
end

function get_export_fieldnames(o::Dict{Symbol,Any})
    local export_fieldnames::NTuple{N,Symbol} where {N}
    if haskey(o, :write)
        export_fieldnames = o[:write]
    else
        export_fieldnames = DEFAULT_EXPORT_FIELDS
    end
    return export_fieldnames
end

function export_results(dh::AbstractDataHandler, options::ExportOptions, chunk_id::Int,
                        timestep::Int, time::Float64)
    options.exportflag || return nothing
    if mod(timestep, options.freq) == 0
        _export_results(dh.chunks[chunk_id], chunk_id, dh.n_chunks, options, timestep, time)
    end
    return nothing
end

function export_reference_results(dh::AbstractDataHandler, options::ExportOptions)
    options.exportflag || return nothing
    @threads :static for chunk_id in eachindex(dh.chunks)
        _export_results(dh.chunks[chunk_id], chunk_id, dh.n_chunks, options, 0, 0.0)
    end
    return nothing
end

function _export_results(b::AbstractBodyChunk, chunk_id::Int, n_chunks::Int,
                         options::ExportOptions, n::Int, t::Float64)
    filename = joinpath(options.vtk, @sprintf("timestep_%05d", n))
    position = get_loc_position(b)
    pvtk_grid(filename, position, b.cells; part=chunk_id, nparts=n_chunks) do vtk
        for (field, type) in options.fields
            vtk[string(field), VTKPointData()] = get_export_field(b.store, field, type)
        end
        vtk["time", VTKFieldData()] = t
    end
    return nothing
end

@inline function get_export_field(s::AbstractStorage, name::Symbol, ::Type{V}) where {V}
    export_field::V = getfield(s, name)
    return export_field
end

@inline function get_loc_position(b::AbstractBodyChunk)
    return @views b.store.position[:, 1:b.ch.n_loc_points]
end
