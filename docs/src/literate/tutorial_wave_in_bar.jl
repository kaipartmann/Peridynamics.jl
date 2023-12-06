# # Pressure wave through bar

# First, import the Peridynamics.jl package:
using Peridynamics

# To start off, we define the geometry of the simulated bar:
# - Length $l_x = 0.2\,\text{m}$
# - Width $l_y = 0.05\,\text{m}$
# - Height $l_z = 0.05\,\text{m}$
# With the point spacing $Δx = \frac{l_x}{60}$, we can create a point cloud.
length_x = 0.2 # [m]
length_y = 0.05 # [m]
length_z = 0.05 # [m]
Δx = length_x / 60 # point spacing
pc = PointCloud(length_x, length_y, length_z, Δx)

# Now we characterize the material model as bond-based with the 
# - Horizon $ δ = 3.015 \cdot Δx$ 
# and the material parameters
# - Density $ρ = 7850\,\mathrm{kg}\,\mathrm{m}^{-3}$
# - Young's modulus $E = 2  10 \, \mathrm{GPa}$
# - Griffith's parameter $G_c = 1000 \, \mathrm{N} \, \mathrm{m}^{-1}$ 
mat = BondBasedMaterial(;
    horizon=3.015Δx, # [m]
    rho=7850.0, # [kg/m^3]
    E=210e9, # [N/m^2]
    Gc=1000.0, # [N/m]
)

# Then we define the boundary conditions, which apply a velocity to the body.
# Since the velocity is supposed to be applied on the left end of the bar, we first
# define the set of points that represent the outermost layer in x-direction of the
# points contained in the point cloud $pc$.
point_set_bc =  findall(pc.position[1,:] .< - length_x / 2 + 1.2 * Δx)
# The applied velocity is represented by one sine oscillation with the period
# $T = 0.00002\,\mathrm{s}$ and the amplitude 
# $v_\mathrm{max} = 0.5\,\mathrm{m}\,\mathrm{s}^{-1}$.
# Afterwards, the velocity is $0$. (see Figure)
# ![](../assets/impact_velocity.png)
T = 0.00002 # [s]
vmax = 0.5 # [m/s]
vwave(t) = t < T ? vmax * sin(2π/T * t) : 0
# Now we submit these conditions. The number indicates the dimension in which the velocity
# is applied, which in this case is the x-dimension signified by 1.
boundary_conditions = [VelocityBC(vwave, point_set_bc, 1)]

# Then we define the time discretization. We calculate 2000 time steps with the
# Velocity Verlet algorithm:
td = VelocityVerlet(2000)

# We assign a name and specify the export settings such as the folder that the data
# will be saved in. Here, 10 determines that every 10th time step
# will be exported.
simulation_name = "PressureWave"
resfolder = joinpath(@__DIR__, "results", simulation_name)
mkpath(resfolder)
es = ExportSettings(resfolder, 10)

# We define a job and pass along the previously set attributes.
job = PDSingleBodyAnalysis(;
    name=simulation_name,
    pc=pc,
    mat=mat,
    bcs=boundary_conditions,
    td=td,
    es=es,
)

# The last step is submitting the job to start the simulation.
#md # ```julia
#md # results = submit(job);
#md # ```

