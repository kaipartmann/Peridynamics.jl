function Base.show(io::IO, ::MIME"text/plain", pdsba::T) where {T<:AbstractPDAnalysis}
    print(io, typeof(pdsba), ": ", pdsba.name)
    return nothing
end

@doc raw"""
    PDSingleBodyAnalysis{T<:AbstractPDMaterial} <: AbstractPDAnalysis

Peridynamic single body analysis.

# Fields
- `name::String`: simulation name
- `pc::PointCloud`: point cloud
- `mat::T`: material model for the body
- `precracks::Vector{PreCrack}`: predefined cracks
- `bcs::Vector{<:AbstractBC}`: boundary conditions
- `ics::Vector{<:AbstractIC}`: initial conditions
- `td::TimeDiscretization`: time discretization
- `es::ExportSettings`: export settings
"""
struct PDSingleBodyAnalysis{T<:AbstractPDMaterial} <: AbstractPDAnalysis
    name::String
    pc::PointCloud
    mat::Union{T,MultiMaterial{T}}
    precracks::Vector{PreCrack}
    bcs::Vector{<:AbstractBC}
    ics::Vector{<:AbstractIC}
    td::TimeDiscretization
    es::ExportSettings
    function PDSingleBodyAnalysis(;
        name::String,
        pc::PointCloud,
        mat::Union{AbstractPDMaterial,MultiMaterial},
        precracks::Vector{PreCrack}=Vector{PreCrack}(),
        bcs::Vector{<:AbstractBC}=Vector{AbstractBC}(),
        ics::Vector{<:AbstractIC}=Vector{AbstractIC}(),
        td::TimeDiscretization,
        es::ExportSettings,
    )
        es.resultfile_prefix = joinpath(es.path, name)
        es.logfile = es.resultfile_prefix * ".log"
        if !es.exportflag
            es.exportfreq = td.n_timesteps + 1
        end
        return new{eltype(mat)}(name, pc, mat, precracks, bcs, ics, td, es)
    end
end

"""
    submit(sim::T) where {T<:AbstractPDAnalysis}

Submit a simulation job to determine the results of the specified analysis.

# Arguments
- `sim::T where {T<:AbstractPDAnalysis}`: simulation job
"""
function submit end

function submit(sim::PDSingleBodyAnalysis)
    timingsummary = @timed begin
        log_header(sim.es)
        log_describe_sim(sim.name, sim.pc, sim.mat, sim.es)
        body = create_simmodel(sim.mat, sim.pc)
        log_describe_interactions(body, sim.es)
        for precrack in sim.precracks
            define_precrack!(body, precrack)
        end
        update_thread_cache!(body)
        calc_damage!(body)
        if sim.td.Δt < 0.0 && sim.td.alg !== :dynrelax
            sim.td.Δt = calc_stable_timestep(body, sim.mat)
        end
        log_simsetup(sim)
        apply_ics!(body, sim.ics)
        if sim.es.exportflag
            export_vtk(body, sim.es.resultfile_prefix, 0, 0.0)
        end
        if sim.td.alg == :verlet
            velocity_verlet!(body, sim)
        elseif sim.td.alg == :dynrelax
            dynamic_relaxation_finite_difference!(body, sim)
        end
    end
    log_closesim(body, timingsummary, sim.es)
    return body
end
