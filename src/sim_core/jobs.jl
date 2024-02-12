
struct Job{S<:SpatialSetup,T<:AbstractTimeSolver}
    spatial_setup::S
    time_solver::T
    es::ExportSettings
    function Job(spatial_setup::S, time_solver::T; path::AbstractString="",
                 writefields::Vector{Symbol}=WRITEFIELDS_DEFAULTS,
                 freq::Int=0) where {S<:SpatialSetup,T<:AbstractTimeSolver}
        es = ExportSettings(spatial_setup, path, freq, writefields)
        new{S,T}(spatial_setup, time_solver, es)
    end
end
