function Base.show(io::IO, ::MIME"text/plain", pdsba::T) where {T <: AbstractPDAnalysis}
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
struct PDSingleBodyAnalysis{M <: AbstractPDMaterial, T <: AbstractTimeDiscretization} <:
       AbstractPDAnalysis
    name::String
    pc::PointCloud
    mat::Union{M, MultiMaterial{M}}
    precracks::Vector{PreCrack}
    bcs::Vector{<:AbstractBC}
    ics::Vector{<:AbstractIC}
    td::T
    es::ExportSettings
    function PDSingleBodyAnalysis(; name::String, pc::PointCloud, mat::PDMaterial{M},
                                  precracks::Vector{PreCrack} = Vector{PreCrack}(),
                                  bcs::Vector{<:AbstractBC} = Vector{AbstractBC}(),
                                  ics::Vector{<:AbstractIC} = Vector{AbstractIC}(), td::T,
                                  es::ExportSettings) where {M <: AbstractPDMaterial,
                                                             T <:
                                                             AbstractTimeDiscretization}
        es.resultfile_prefix = joinpath(es.path, name)
        es.logfile = es.resultfile_prefix * ".log"
        if !es.exportflag
            es.exportfreq = typemax(Int)
        end
        return new{M, T}(name, pc, mat, precracks, bcs, ics, td, es)
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
        init_time_discretization!(sim.td, body, sim.mat)
        log_simsetup(sim)
        time_loop!(body, sim.td, sim.mat, sim.bcs, sim.ics, sim.es)
    end
    log_closesim(body, timingsummary, sim.es)
    return body
end
