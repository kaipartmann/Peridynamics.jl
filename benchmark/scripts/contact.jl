using Peridynamics

const RESPATH = length(ARGS) ≥ 1 ? ARGS[1] : "results"
const SIMNAME = length(ARGS) ≥ 2 ? ARGS[2] : "contact"
const EXPORT_VTK = length(ARGS) ≥ 3 ? parse(Bool, ARGS[3]) : true

function cube_impact(lp, Δxp, lc, Δxc, v0, nt)
    pcp = PointCloud(lp, lp, 0.05lp, Δxp)
    matp = BondBasedMaterial(horizon=3.015Δxp, rho=2e-6, E=2e4, Gc=0.01)
    bsplate = BodySetup(pcp, matp)
    pcc = PointCloud(lc, lc, lc, Δxc)
    pcc.position[3,:] .+= lc/2 + 0.025lp + 0.01Δxc
    matc = BondBasedMaterial(horizon=3.015Δxc, rho=8e-6, E=2.1e5, Gc=3)
    ics = [VelocityIC(v0, 1:pcc.n_points, 3)]
    bscube = BodySetup(pcc, matc; ics=ics)
    bodies = [bsplate, bscube]
    c = [Contact((1,2), max(Δxp, Δxc), 1e6)]
    vv = VelocityVerlet(nt)
    path = joinpath(RESPATH, SIMNAME)
    ispath(path) && rm(path; recursive=true, force=true)
    !ispath(path) && mkpath(path) # create the path if it does not exist
    es = EXPORT_VTK ? ExportSettings(path, 10) : ExportSettings(path, typemax(Int))
    job = PDContactAnalysis(name=SIMNAME, body_setup=bodies, contact=c, td=vv, es=es)
    submit(job)
    return nothing
end

println()
println("-"^70)
@show SIMNAME
@show EXPORT_VTK
println("-"^70)
println()

# @time cube_impact(lp, Δxp, lc, Δxc, v0, nt)
@time cube_impact(1, 1/30, 0.1, 0.1/5, -1000, 1000)
@time cube_impact(1, 1/150, 0.1, 0.1/15, -1000, 2000)
