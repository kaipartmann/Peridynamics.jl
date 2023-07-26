# # [Contact analysis](@id howto_contact_analysis)
# Short guide on how to run a peridynamic analysis with contact of multiple bodies.

# ## Steps
# ### 1. Define the bodies
# To define the bodies, create [`BodySetup`](@ref) instances. Each `BodySetup` is created
# with [step 1-4 from single body analysis guide](@ref howto_single_body_analysis).
# Optionally choose if a body should contribute to the automatic calculation of the stable
# time step with the `calc_timestep` keyword. Be careful if you use this option!
# ### 2. Define contact
# Define the contact of two bodies by creating a [`Contact`](@ref) instance.
# ### 3. Choose temporal discretization settings
# Create a [`VelocityVerlet`](@ref) instance.
# ### 4. Choose export settings
# Create a [`ExportSettings`](@ref) instance.
# ### 5. Create the simulation job
# Create a [`PDContactAnalysis`](@ref) instance and choose a simulation name.
# ### 6. Run the simulation job
# Run the simulation with the [`submit`](@ref) function.

# ## Examples
# Please refer to the [Peridynamics.jl logo tutorial](@ref tutorial_logo).
