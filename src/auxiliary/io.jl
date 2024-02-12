
struct ExportSettings
    exportflag::Bool
    root_path::String
    results_path::String
    logfile_path::String
    exportfreq::Int
    fields::Vector{Symbol}

    function ExportSettings(spatial_setup::SpatialSetup, path::AbstractString, freq::Int,
                            writefields::Vector{Symbol})
        if isempty(path) && freq == 0
            #TODO: Warning if writefields != defaults
            return new(false, "", "", "", 0, Symbol[])
        else
            results_path = joinpath(path, "vtk")
            logfile_path = joinpath(path, "logfile.log")
            return new(true, path, results_path, logfile_path, freq, writefields)
        end
    end
end
