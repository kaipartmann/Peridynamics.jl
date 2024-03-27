const EXPORT_KWARGS = (:path, :freq, :fields)
const DEFAULT_EXPORT_FIELDS = (:displacement, :damage)

struct ExportOptions
    exportflag::Bool
    root::String
    vtk::String
    logfile::String
    freq::Int
    fields::Vector{Symbol}
end

function ExportOptions(root::String, freq::Int, fields::Vector{Symbol})
    vtk = joinpath(root, "vtk")
    logfile = joinpath(root, "logfile.log")
    return ExportOptions(true, root, vtk, logfile, freq, fields)
end

function ExportOptions()
    return ExportOptions(false, "", "", "", 0, Vector{Symbol}())
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

    fields = get_export_fields(S, o)

    if isempty(root)
        eo = ExportOptions()
    else
        eo = ExportOptions(root, freq, fields)
    end

    return eo
end

function get_export_fields(::Type{S}, o::Dict{Symbol,Any}) where {S}
    local _fields::NTuple{N,Symbol} where {N}
    if haskey(o, :fields)
        _fields = o[:fields]
    else
        _fields = DEFAULT_EXPORT_FIELDS
    end

    fields = [field for field in _fields]
    check_storage_fields(S, fields)

    return fields
end

function check_storage_fields(::Type{S}, fields::Vector{Symbol}) where {S}
    allowed_fields = fieldnames(S)
    for f in fields
        if !in(f, allowed_fields)
            msg = "unknown field $(name) specified for export!\n"
            msg *= "Allowed fields for $S:\n"
            for allowed_name in fieldnames(S)
                msg *= "  - $allowed_name\n"
            end
            throw(ArgumentError(msg))
        end
    end
    return nothing
end

function _export_results(b::AbstractBodyChunk, chunk_id::Int, n_chunks::Int,
                         options::ExportOptions, n::Int, t::Float64)
    filename = @sprintf("timestep_%06d", n)
    position = get_loc_position(b)
    pvtk_grid(filename, position, b.cells; part=chunk_id, nparts=n_chunks) do vtk
        for field in options.fields
            vtk[string(field), VTKPointData()] = get_storage_field(b.store, field)
        end
        vtk["time", VTKFieldData()] = t
    end
    return nothing
end

@inline function get_loc_position(b::AbstractBodyChunk)
    return @views b.store.position[:, 1:b.ch.n_loc_points]
end

function init_export(options::ExportOptions)
    options.exportflag || return ""
    run_pwd = pwd()
    mkpath(options.vtk)
    cd(options.vtk)
    return run_pwd
end

function finish_export(options::ExportOptions, run_pwd::AbstractString)
    options.exportflag && cd(run_pwd)
    return nothing
end
