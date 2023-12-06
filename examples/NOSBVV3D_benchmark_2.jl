using Peridynamics
using LinearAlgebra
# BLAS.set_num_threads(1)

RESPATH::String = length(ARGS) ≥ 1 ? ARGS[1] : joinpath("results", "NOSBVV3D")

function benchmark_2_xwave(simname, l, lyz, Δx, vmax, T, nt)
    pc = PointCloud(l, lyz, lyz, Δx)
    δ = 3.015Δx
    mat = NOSBMaterial(horizon=δ,rho=8000,E=210e9,nu=0.25,epsilon_c=0.01,stabilization=0)
    points_fix_left =  findall(pc.position[1,:] .< - l / 2 + 1.2 * Δx)
    vwave(t) = t < T ? vmax * sin(2π/T * t) : 0
    bcs = [VelocityBC(vwave, points_fix_left, 1)]
    td = VelocityVerlet(nt)
    path = joinpath(RESPATH, simname)
    ispath(path) && rm(path; recursive=true, force=true)
    !ispath(path) && mkpath(path) # create the path if it does not exist
    es = ExportSettings(path, 10)
    job = PDSingleBodyAnalysis(name=simname, pc=pc, mat=mat, bcs=bcs, td=td, es=es)
    submit(job)
    return nothing
end

##-- SIMULATION PARAMETERS
simname = "NOSBVV3D_benchmark_2_xwave"
l = 0.2
lyz = 0.01l
npyz = 4
Δx = lyz / npyz
vmax = 10
Δt = 7.05e-7
T = 30Δt
nt = 2000

##-- PRECOMPILATION
_lyz = 0.5l
_npyz = 5
_Δx = _lyz / _npyz
_nt = 500
printstyled("----- PRECOMPILATION -----\n"; bold=true, color=:red)
@time benchmark_2_xwave(simname, l, _lyz, _Δx, vmax, T, _nt)

##-- SIMULATION
printstyled("----- SIMULATION -----\n"; bold=true, color=:red)
@time benchmark_2_xwave(simname, l, lyz, Δx, vmax, T, nt)
