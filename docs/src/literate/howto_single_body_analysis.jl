# # [Single body analysis](@id howto_single_body_analysis)
# Short guide on how to run a peridynamics analysis with only one body.

# ## Steps
# ### 1. Define a point cloud
# Please refer to the [guide on how to use point clouds](@ref howto_pointclouds).
# ### 2. Set material properties
# Please refer to the
# [guide on how to use material formulations](@ref howto_matformulations).
# ### 3. Create predefined cracks [optional]
# Please refer to the [guide on how to use predefined cracks](@ref howto_precracks).
# ### 4. Assign boundary and initial conditions [optional]
# ##### Boundary conditions, updated every time step
# - Velocity boundary condition, [`VelocityBC`](@ref)
# - Force density boundary condition [`ForceDensityBC`](@ref)
# - Position dependend velocity boundary condition, [`PosDepVelBC`](@ref)
# ##### Initial conditions, only assigned before time-loop
# - Velocity initial condition [`VelocityIC`](@ref)
# ### 5. Assign a temporal discretization
# Create either a [`VelocityVerlet`](@ref) or a [`DynamicRelaxation`](@ref) instance.
# ### 6. Choose export settings
# Create a [`ExportSettings`](@ref) instance.
# ### 7. Create the simulation job
# Create a [`PDSingleBodyAnalysis`](@ref) instance and choose a simulation name.
# ### 8. Run the simulation job
# Run the simulation with the [`submit`](@ref) function.

# ## Examples
# Please refer to the following tutorials:
#   - [Mode I tension quasi-static](@ref tutorial_tension_static)
#   - [Mode I tension dynamic](@ref tutorial_tension_dynfrac)
#   - [Mode I tension dynamic with predefined crack](@ref tutorial_tension_precrack)
