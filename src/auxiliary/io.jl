const ExportField = Tuple{Symbol,DataType}

struct ExportOptions{F<:Function}
    exportflag::Bool
    root::String
    vtk::String
    logfile::String
    freq::Int
    eff::F

    function ExportOptions(root::String, freq::Int, eff::F) where {F<:Function}
        if isempty(root)
            _eff = ()->()
            return new{typeof(_eff)}(false, "", "", "", 0, _eff)
        end
        vtk = joinpath(root, "vtk")
        logfile = joinpath(root, "logfile.log")
        return new{F}(true, root, vtk, logfile, freq, eff)
    end
end

function (eo::ExportOptions{F})(s::AbstractStorage) where {F<:Function}
    return eo.eff(s)
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

    eff = get_exported_fields_function(S, o)

    return ExportOptions(root, freq, eff)
end

function export_disp_and_dmg(s::AbstractStorage)
    return (("displacement", s.displacement), ("damage", s.damage))
end

macro eff(vars...)
    local args = Expr[]
    for var in vars
        var isa QuoteNode || error("only symbols are allowed args!\n")
        local name = string(var)[begin+1:end]
        push!(args, :($(name), Core.getfield(x, $(esc(var)))))
    end
    return Expr(:->, :x, Expr(:tuple, args...))
end

function get_exported_fields_function(::Type{S}, o::Dict{Symbol,Any}) where {S}
    if haskey(o, :eff)
        eff = o[:eff]
    else
        eff = export_disp_and_dmg
    end
    return eff
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
        for (name, field) in options(b.store)
            vtk[name] = field
        end
        vtk["time", VTKFieldData()] = t
    end
    return nothing
end

@inline function get_loc_position(b::AbstractBodyChunk)
    return @views b.store.position[:, 1:b.ch.n_loc_points]
end
