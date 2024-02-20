
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

function get_export_options(::Type{M}, o::Dict{Symbol,Any}) where {M<:AbstractMaterial}
    local root::String
    local freq::Int
    local fields::NTuple{N,Symbol} where {N}

    if haskey(o, :path)
        root = string(o[:path])
    else
        root = ""
    end

    if haskey(o, :freq)
        freq = Int(o[:freq])
    else
        freq = 0
    end
    freq < 0 && throw(ArgumentError("`freq` should be larger than zero!\n"))

    if haskey(o, :write)
        fields = o[:write]
    else
        fields = default_export_fields(M)
    end

    return ExportOptions(root, freq, fields)
end

function export_results(dh::AbstractDataHandler, options::ExportOptions, timestep::Int,
                        time::Float64)
    if mod(timestep, options.freq) == 0
        _export_results(dh, options, timestep, time)
    end
    return nothing
end
