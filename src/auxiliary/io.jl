struct ExportOptions{N}
    exportflag::Bool
    root::String
    vtk::String
    logfile::String
    freq::Int
    fields::NTuple{N,Symbol}

    function ExportOptions(root::String, freq::Int, fields::NTuple{N,Symbol}) where {N}
        if isempty(root)
            return new{0}(false, "", "", "", 0, NTuple{0,Symbol}())
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
        root = abspath(string(o[:path]))
        freq = Int(o[:freq])
    elseif haskey(o, :path) && !haskey(o, :freq)
        root = abspath(string(o[:path]))
        freq = 10
    elseif !haskey(o, :path) && haskey(o, :freq)
        msg = "if `freq` is spedified, the keyword `path` is also needed!\n"
        throw(ArgumentError(msg))
    else
        root = ""
        freq = 0
    end
    freq < 0 && throw(ArgumentError("`freq` should be larger than zero!\n"))

    exfields = get_exported_fields(S, o)

    return ExportOptions(root, freq, exfields)
end

function get_exported_fields(::Type{S}, o::Dict{Symbol,Any}) where {S}
    if haskey(o, :fields)
        fields = o[:fields]
    else
        fields = DEFAULT_EXPORT_FIELDS
    end
    return fields
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
    filename = @sprintf("timestep_%05d", n)
    position = get_loc_position(b)
    pvtk_grid(filename, position, b.cells; part=chunk_id, nparts=n_chunks) do vtk
        for field in options.fields
            vtk[string(field)] = get_storage_field(b.store, field)
        end
        vtk["time", VTKFieldData()] = t
    end
    return nothing
end

@inline function get_loc_position(b::AbstractBodyChunk)
    return @views b.store.position[:, 1:b.ch.n_loc_points]
end
