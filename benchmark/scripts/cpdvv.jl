using Peridynamics

const RESPATH = length(ARGS) ≥ 1 ? ARGS[1] : "results"
const SIMNAME = length(ARGS) ≥ 2 ? ARGS[2] : "cpdvv"
const EXPORT_VTK = length(ARGS) ≥ 3 ? parse(Bool, ARGS[3]) : true

function mode_I_tension_dynamic(lx, ly, lz, Δx, v0, nt)
    pc = PointCloud(lx, ly, lz, Δx)
    δ = 3.015Δx
    mat = ContinuumBasedMaterial(horizon=δ, rho=8e-6, E=2.1e5, nu=0.3, Gc=0.05)
    a = 0.5lx
    set_a = findall(p -> p[1] ≤ -lx/2+a && 0 ≤ p[2] ≤ 2δ, eachcol(pc.position))
    set_b = findall(p -> p[1] ≤ -lx/2+a && -2δ ≤ p[2] < 0, eachcol(pc.position))
    precrack = PreCrack(set_a, set_b)
    set_top = findall(p -> p[2] ≥ ly/2-Δx, eachcol(pc.position))
    set_bottom = findall(p -> p[2] ≤ -ly/2+Δx, eachcol(pc.position))
    bc_bottom = VelocityBC(t -> -v0, set_bottom, 2)
    bc_top = VelocityBC(t -> v0, set_top, 2)
    vv = VelocityVerlet(nt)
    path = joinpath(RESPATH, SIMNAME)
    ispath(path) && rm(path; recursive=true, force=true)
    !ispath(path) && mkpath(path) # create the path if it does not exist
    es = EXPORT_VTK ? ExportSettings(path, 10) : ExportSettings(path, typemax(Int))
    job = PDSingleBodyAnalysis(name=SIMNAME, pc=pc, mat=mat, precracks=[precrack],
                               bcs=[bc_bottom,bc_top], td=vv, es=es)
    submit(job)
    return nothing
end

println()
println("-"^70)
@show SIMNAME
@show EXPORT_VTK
println("-"^70)
println()

@time mode_I_tension_dynamic(1.0, 1.0, 1.0, 1/5, 10, 1000) # compilation run
@time mode_I_tension_dynamic(1.0, 1.0, 0.1, 1/40, 10, 1000)
