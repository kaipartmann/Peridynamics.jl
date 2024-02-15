
struct Job{S<:SpatialSetup,T<:AbstractTimeSolver}
    spatial_setup::S
    time_solver::T
    options::ExportOptions

    function Job(spatial_setup::S, time_solver::T, options::ExportOptions) where {S,T}
        pre_submission_check(spatial_setup)
        return new{S,T}(spatial_setup, time_solver, options)
    end
end

function Job(spatial_setup::S, time_solver::T; kwargs...) where {S,T}
    o = Dict{Symbol,Any}(kwargs)
    check_kwargs(o, JOB_KWARGS)
    options = get_export_options(material_type(spatial_setup), o)
    Job(spatial_setup, time_solver, options)
end
