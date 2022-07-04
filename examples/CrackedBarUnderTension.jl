#===========================================================================================
    BAR UNDER TENSION - EXAMPLE SCRIPT
===========================================================================================#

##--   1. LOAD MODULES   -------------------------------------------------------------------
using Peridynamics

##--   2. DEFINE THE POINTCLOUD   ----------------------------------------------------------
length_x = 0.05 # [m]
length_y = 0.05 # [m]
length_z = 0.005 # [m]
Δx = length_y / 60 # point spacing: 60 points along y-describe_interactions
pc = PointCloud(length_x, length_y, length_z, Δx)

##--   3. DEFINE THE MATERIAL   ------------------------------------------------------------
mat = BondBasedMaterial(;
    horizon = 3.015Δx, # [m]
    rho = 7850.0, # [kg/m^3]
    E = 210e9, # [N/m^2]
    Gc = 1000.0, # [N/m]
)

##--   3. PREDEFINED CRACKS   --------------------------------------------------------------
cracklength = 0.5 * length_x
precrack_set_a = findall((pc.position[2, :] .>= 0) .& (pc.position[2, :] .< 12 * Δx) .&
                         (pc.position[1, :] .<= -length_x / 2 + cracklength))
precrack_set_b = findall((pc.position[2, :] .<= 0) .& (pc.position[2, :] .> -12 * Δx) .&
                         (pc.position[1, :] .<= -length_x / 2 + cracklength))
precracks = [PreCrack(precrack_set_a, precrack_set_b)]

##--   4. DEFINE INITIAL AND BOUNDARY CONDITIONS   -----------------------------------------
bc_set_top = findall(pc.position[2,:] .> length_y / 2 - 5.1 * Δx)
bc_set_bottom = findall(pc.position[2,:] .< -length_y / 2 + 5.1 * Δx)
bc_top = VelocityBC(t -> 0.1, bc_set_top, 2)
bc_bottom = VelocityBC(t -> -0.1, bc_set_bottom, 2)
boundary_conditions = [bc_top, bc_bottom]

##--   5. DEFINE TIME DISCRETIZATION   -----------------------------------------------------
td = TimeDiscretization(2000)

##--   6. DEFINE EXPORT SETTINGS   ---------------------------------------------------------
simulation_name = "CrackedBarUnderTension"
resfolder = joinpath(@__DIR__, "results", simulation_name)
mkpath(resfolder)
es = ExportSettings(resfolder, 10)

##--   7. DEFINE THE JOB   -----------------------------------------------------------------
job = PDSingleBodyAnalysis(
    name = simulation_name,
    pc = pc,
    mat = mat,
    precracks = precracks,
    bcs = boundary_conditions,
    td = td,
    es = es,
)

##--   8. SUBMIT THE JOB / RUN THE ANALYSIS   ----------------------------------------------
results = submit(job);
