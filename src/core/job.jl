
struct Job{S<:SpatialSetup,T<:AbstractTimeSolver}
    spatial_setup::S
    time_solver::T
    opt::ExportOptions
end

function Job(spatial_setup::S, time_solver::T; kwargs...) where {S,T}
    o = Dict{Symbol,Any}(kwargs)
    check_kwargs(o, JOB_KWARGS)
    opt = get_export_options(material_type(spatial_setup), o)
    Job(spatial_setup, time_solver, opt)
end
