#===========================================================================================
    IMPACT-SIMULATION FOR THE LOGO OF THE PERIDYNAMICS PACKAGE
===========================================================================================#

##--   1. LOAD MODULES   -------------------------------------------------------------------
using Peridynamics

##--   2. DEFINE THE POINTCLOUDS   ---------------------------------------------------------
# PLATE - INDEX: p
dimXYₚ = 0.1 # [m]
dimZₚ = 0.01 # [m]
Δxₚ = dimXYₚ / 50 # [m]
pcₚ = PointCloud(dimXYₚ, dimXYₚ, dimZₚ, Δxₚ)

# SPHERES OF LOGO - INDEX: s
Øₛ = 0.03 # [m]
Vₛ = 4 / 3 * π * (Øₛ / 2)^3 # [m³]
Δxₛ = Øₛ / 20 # [m]
pcₛ₀ = PointCloud(Øₛ, Øₛ, Øₛ, Δxₛ)
sphere_point_set = @views findall(
    sqrt.(
        pcₛ₀.position[1, :] .^ 2 + pcₛ₀.position[2, :] .^ 2 + pcₛ₀.position[3, :] .^ 2
    ) .<= Øₛ / 2,
)
pcₛ₀.position[3, sphere_point_set] .+= Øₛ / 2 + dimZₚ / 2 + 1.1Δxₚ
pcₛ₁ = PointCloud(
    pcₛ₀.position[:, sphere_point_set],
    pcₛ₀.volume[sphere_point_set],
    zeros(Bool, length(sphere_point_set)),
    pcₛ₀.radius[sphere_point_set],
    length(sphere_point_set),
)
pcₛ₂ = deepcopy(pcₛ₁)
pcₛ₃ = deepcopy(pcₛ₁)
rₗ = Øₛ / 2 + 0.2 * Øₛ
pcₛ₁.position[2, :] .+= rₗ
pcₛ₂.position[1, :] .+= rₗ * cos(30π / 180)
pcₛ₂.position[2, :] .-= rₗ * sin(30π / 180)
pcₛ₃.position[1, :] .-= rₗ * cos(30π / 180)
pcₛ₃.position[2, :] .-= rₗ * sin(30π / 180)

##--   3. DEFINE THE MATERIAL   ------------------------------------------------------------
matₚ = BondBasedMaterial(;
    point_spacing=Δxₚ,
    horizon=3.015Δxₚ,
    density=2000.0,
    young_modulus=30e9,
    critical_energy_release_rate=10.0,
)
matₛ = BondBasedMaterial(;
    point_spacing=Δxₛ,
    horizon=3.015Δxₛ,
    density=7850.0,
    young_modulus=210e9,
    critical_energy_release_rate=1000.0,
)

##--   4. DEFINE INITIAL CONDITIONS   ------------------------------------------------------
icₛ = [VelocityIC(-20.0, 1:(pcₛ₁.n_points), 3)]

##--   5. DEFINE THE BODY-SETUP   ----------------------------------------------------------
plate = BodySetup(pcₚ, matₚ)
sphere1 = BodySetup(pcₛ₁, matₛ; ics=icₛ)
sphere2 = BodySetup(pcₛ₂, matₛ; ics=icₛ, calc_timestep=false)
sphere3 = BodySetup(pcₛ₃, matₛ; ics=icₛ, calc_timestep=false)
body_setup = [plate, sphere1, sphere2, sphere3]

##--   6. DEFINE CONTACT   -----------------------------------------------------------------
contact = [Contact((1, 2), Δxₚ), Contact((1, 3), Δxₚ), Contact((1, 4), Δxₚ)]

##--   7. DEFINE TIME DISCRETIZATION   -----------------------------------------------------
td = TimeDiscretization(3000)

##--   8. DEFINE EXPORT SETTINGS   ---------------------------------------------------------
simulation_name = "Logo"
resfolder = joinpath(@__DIR__, "results", simulation_name)
mkpath(resfolder)
es = ExportSettings(resfolder, 10)

##--   9. DEFINE THE JOB   -----------------------------------------------------------------
JOB = PDContactAnalysis(;
    name=simulation_name, body_setup=body_setup, contact=contact, td=td, es=es
)

##--   11. SUBMIT THE JOB FOR ANALYSIS   ---------------------------------------------------
results = submit(JOB);
