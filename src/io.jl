"""
    ExportSettings

Export settings.

# Fields
- `path::String`: path where results will be saved
- `exportfreq::Int`: export frequency, will export every `exportfreq`-th timestep
- `resultfile_prefix::String`: prefix of the result-filename
- `logfile::String`: name of logfile
- `exportflag::Bool`: disable export for a simulation where saved results are not needed

---
```julia
ExportSettings([path::String, freq::Int])
```

Create `ExportSettings` only by `path` and `freq`. If no arguments are specified, the
`exportflag` will be set to `false` and export disabled.

# Arguments
- `path::String`: path where results will be saved
- `freq::Int`: export frequency
"""
mutable struct ExportSettings
    path::String
    exportfreq::Int
    resultfile_prefix::String
    logfile::String
    exportflag::Bool
end

ExportSettings() = ExportSettings("", 0, "", "", false)
ExportSettings(path::String, freq::Int) = ExportSettings(path, freq, "", "", true)

function Base.show(io::IO, ::MIME"text/plain", es::ExportSettings)
    print(io, typeof(es))
    return nothing
end


function log_header(es::ExportSettings)
    msg = "="^70 * "\n"
    msg *= string("PERIDYNAMIC SIMULATION ON ", nthreads(), " THREADS\n")
    msg *= "="^70 * "\n"
    print(msg)
    if es.exportflag
        open(es.logfile, "w") do io
            return write(io, msg)
        end
    end
    return nothing
end

function log_describe_sim(
    simname::String, pc::PointCloud, mat::AbstractPDMaterial, es::ExportSettings
)
    msg = "Peridynamic single body analysis: " * simname * "\nMaterial parameters:\n"
    msg *= @sprintf("  - Number of material points [-]:      %30g\n", pc.n_points)
    msg *= describe_mat(mat)
    print(msg)
    if es.exportflag
        open(es.logfile, "a") do io
            return write(io, msg)
        end
    end
    return nothing
end

function describe_mat(mat::AbstractPDMaterial)
    msg = @sprintf "  - Material model: %50s\n" typeof(mat)
    msg *= @sprintf "  - Horizon δ [m]:                      %30g\n" mat.δ
    msg *= @sprintf "  - Density ρ [kg/m³]:                  %30g\n" mat.rho
    msg *= @sprintf "  - Young's modulus E [N/m²]:           %30g\n" mat.E
    msg *= @sprintf "  - Poisson ratio ν [-]:                %30g\n" mat.nu
    msg *= @sprintf "  - Critical bond stretch εc[-]:        %30g\n" mat.εc
    return msg
end

function log_describe_interactions(part::AbstractPDBody, es::ExportSettings)
    msg = "Interactions:\n"
    msg *= describe_interactions(part)
    msg *= @sprintf(
        "Total memory used by body [MB]:         %30g\n", Base.summarysize(part) * 1e-6
    )
    print(msg)
    if es.exportflag
        open(es.logfile, "a") do io
            return write(io, msg)
        end
    end
    return nothing
end

function describe_interactions(part::AbstractPDBody)
    msg = @sprintf "  - Number of bonds [-]:                %30d\n" part.n_bonds
    return msg
end

function log_simsetup(sim::AbstractPDAnalysis)
    msg = ""
    if sim.es.exportflag
        msg *= "Export setup:\n"
        msg *= @sprintf "  - Export frequency:                   %30g\n" sim.es.exportfreq
        resfile_basename = basename(sim.es.resultfile_prefix)
        if length(resfile_basename) <= 48
            msg *= @sprintf "  - Export file name: %48s\n" resfile_basename
        else
            msg *= @sprintf "  - Export file name:\n    %s\n" resfile_basename
        end
    end
    timedisc = ""
    if sim.td.alg == :verlet
        timedisc = "dynamic (Velocity-Verlet algorithm)"
    elseif sim.td.alg == :dynrelax
        timedisc = "quasi-static (adaptive dynamic relaxation)"
    end
    msg *= @sprintf "Time discretization: %49s\n" timedisc
    msg *= @sprintf "  - Time step Δt [s]:                   %30g\n" sim.td.Δt
    msg *= @sprintf "  - Number of time steps [-]:           %30g\n" sim.td.n_timesteps
    msg *= @sprintf "  - Simulation time horizon [s]:        %30g\n" sim.td.n_timesteps *
        sim.td.Δt
    print(msg)
    if sim.es.exportflag
        open(sim.es.logfile, "a") do io
            return write(io, msg)
        end
    end
    return nothing
end

function log_closesim(part::AbstractPDBody, timingsummary::NamedTuple, es::ExportSettings)
    msg = @sprintf(
        "✓ Simulation completed after %g seconds\nResults:\n", timingsummary.time
    )
    msg *= @sprintf(
        "  - Max. abs. x-displacement [m]:       %30g\n",
        abs(maximum(part.displacement[1, :]))
    )
    msg *= @sprintf(
        "  - Max. abs. y-displacement [m]:       %30g\n",
        abs(maximum(part.displacement[2, :]))
    )
    msg *= @sprintf(
        "  - Max. abs. z-displacement [m]:       %30g\n",
        abs(maximum(part.displacement[3, :]))
    )
    msg *= @sprintf("  - Max. damage [-]:                    %30g\n", maximum(part.damage))
    print(msg)
    if es.exportflag
        open(es.logfile, "a") do io
            return write(io, msg)
        end
    end
    return nothing
end


function export_vtk(body::AbstractPDBody, expfile::String, timestep::Int, time::Float64)
    filename = expfile * "_t" * string(timestep)
    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (j,)) for j in 1:body.n_points]
    vtkfile = vtk_grid(filename, body.position, cells)
    vtkfile["Damage", VTKPointData()] = body.damage
    vtkfile["ForceDensity", VTKPointData()] = @views body.b_int[:, :, 1]
    vtkfile["Displacement", VTKPointData()] = body.displacement
    vtkfile["Velocity", VTKPointData()] = body.velocity
    vtkfile["Time", VTKFieldData()] = time
    vtk_save(vtkfile)
    return nothing
end
