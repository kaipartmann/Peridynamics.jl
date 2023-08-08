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

function log_describe_sim(simname::String, pc::PointCloud, mat::PDMaterial,
                          es::ExportSettings)
    msg = "Peridynamic single body analysis: " * simname * "\n"
    msg *= "Geometry:\n"
    msg *= describe_geometry(pc)
    msg *= @sprintf("Material Parameters: %49s\n", eltype(mat))
    msg *= describe_mat(mat)
    print(msg)
    if es.exportflag
        open(es.logfile, "a") do io
            return write(io, msg)
        end
    end
    return nothing
end

function describe_geometry(pc::PointCloud)
    minima = minimum(pc.position; dims = 2)
    maxima = maximum(pc.position; dims = 2)
    minmax_x = @sprintf("%.6g / %.6g", minima[1], maxima[1])
    minmax_y = @sprintf("%.6g / %.6g", minima[2], maxima[2])
    minmax_z = @sprintf("%.6g / %.6g", minima[3], maxima[3])
    msg = @sprintf("  - Number of material points [-]:      %30g\n", pc.n_points)
    msg *= @sprintf("  - Min. / Max. x-axis [m]:             %30s\n", minmax_x)
    msg *= @sprintf("  - Min. / Max. y-axis [m]:             %30s\n", minmax_y)
    msg *= @sprintf("  - Min. / Max. z-axis [m]:             %30s\n", minmax_z)
    return msg
end

function describe_mat(mat::AbstractPDMaterial)
    msg = @sprintf "  - Horizon δ [m]:                      %30g\n" mat.δ
    msg *= @sprintf "  - Density ρ [kg/m³]:                  %30g\n" mat.rho
    msg *= @sprintf "  - Young's modulus E [N/m²]:           %30g\n" mat.E
    msg *= @sprintf "  - Poisson ratio ν [-]:                %30g\n" mat.nu
    msg *= @sprintf "  - Griffith's Parameter Gc [N/m]:      %30g\n" mat.Gc
    msg *= @sprintf "  - Critical bond stretch εc[-]:        %30g\n" mat.εc
    return msg
end

function describe_mat(mat::MultiMaterial)
    msg = ""
    for (i, m) in enumerate(mat.materials)
        n_points_mat = length(findall(mat.matofpoint .== i))
        msg *= @sprintf("Parameters for %d points:\n", n_points_mat)
        msg *= describe_mat(m)
    end
    return msg
end

function log_describe_interactions(part::AbstractPDBody, es::ExportSettings)
    mem_used_mb = Base.summarysize(part) * 1e-6
    msg = "Interactions:\n"
    msg *= describe_interactions(part)
    msg *= @sprintf("Total memory used by body [MB]:         %30g\n", mem_used_mb)
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
    msg *= @sprintf "Time discretization: %49s\n" typeof(sim.td)
    msg *= @sprintf "  - Time step Δt [s]:                   %30g\n" sim.td.Δt
    msg *= @sprintf "  - Number of time steps [-]:           %30g\n" sim.td.n_steps
    time_horizon = sim.td.n_steps * sim.td.Δt
    msg *= @sprintf "  - Simulation time horizon [s]:        %30g\n" time_horizon
    print(msg)
    if sim.es.exportflag
        open(sim.es.logfile, "a") do io
            return write(io, msg)
        end
    end
    return nothing
end

function log_closesim(body::AbstractPDBody, timingsummary::NamedTuple, es::ExportSettings)
    msg = @sprintf("✓ Simulation completed after %g seconds\nResults:\n",
                   timingsummary.time)
    msg *= log_displacement_and_damage(body)
    print(msg)
    if es.exportflag
        open(es.logfile, "a") do io
            return write(io, msg)
        end
    end
    return nothing
end

function log_displacement_and_damage(body::AbstractPDBody)
    minima = minimum(body.displacement; dims = 2)
    maxima = maximum(body.displacement; dims = 2)
    minmax_x = @sprintf("%.6g / %.6g", minima[1], maxima[1])
    minmax_y = @sprintf("%.6g / %.6g", minima[2], maxima[2])
    minmax_z = @sprintf("%.6g / %.6g", minima[3], maxima[3])
    msg = ""
    msg *= @views @sprintf("  - Min. / Max. x-displacement [m]:     %30s\n", minmax_x)
    msg *= @views @sprintf("  - Min. / Max. y-displacement [m]:     %30s\n", minmax_y)
    msg *= @views @sprintf("  - Min. / Max. z-displacement [m]:     %30s\n", minmax_z)
    msg *= @sprintf("  - Max. damage [-]:                    %30g\n", maximum(body.damage))
    return msg
end

@timeit TO function export_vtk(body::AbstractPDBody, expfile::String, timestep::Int, time::Float64)
    filename = @sprintf("%s_t%04d", expfile, timestep)
    vtkfile = vtk_grid(filename, body.position, body.cells)
    vtkfile["Damage", VTKPointData()] = body.damage
    vtkfile["ForceDensity", VTKPointData()] = @views body.b_int[:, :, 1]
    vtkfile["Displacement", VTKPointData()] = body.displacement
    vtkfile["Velocity", VTKPointData()] = body.velocity
    vtkfile["Time", VTKFieldData()] = time
    vtk_save(vtkfile)
    return nothing
end

function log_timers(es::ExportSettings)
    if es.exportflag
        open(es.logfile, "a") do io
            write(io, "\n")
            show(IOContext(io, :displaysize => (24,150)), TO)
        end
    else
        println()
        show(TO)
        println()
    end
    return nothing
end
