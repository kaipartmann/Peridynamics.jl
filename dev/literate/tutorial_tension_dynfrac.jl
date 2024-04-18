# # [Mode I tension dynamic with fracture](@id tutorial_tension_dynfrac)

#-
# Import the package:
using Peridynamics
using Peridynamics.AbaqusMeshConverter

#-
# Read and convert the [Abaqus FEM mesh of a tensile test](https://github.com/kaipartmann/Peridynamics.jl/blob/main/docs/src/assets/TensileTestMesh.inp)
# into a point cloud for the peridynamic model:
## insert your correct path to the downloaded mesh file!
pos, vol, pointsets = read_inp(joinpath(@__DIR__, "..", "assets", "TensileTestMesh.inp"))
# ![](../assets/TensileTestMesh.png)

#-
# Create a body with the points from the mesh: (The bond-based material model is used here.)
b = Body(BBMaterial(), pos, vol)

# Convert the element sets defined in Abaqus into point sets of the Body:
point_set!(b, :top, pointsets["top"])
point_set!(b, :bottom, pointsets["bottom"])

# Specify the material parameters as:

# | material parameter | value |
# |:--------|:-------------|
# | Horizon $ δ $ | $0.01 \, \mathrm{m}$ |
# | Density $ρ$ | $ 2700 \,\mathrm{kg}\,\mathrm{m}^{-3}$ |
# | Young's modulus $E$ | $ 70 \cdot 10^{9} \, \mathrm{Pa}$ |
# | Griffith's parameter $G_c$ | $100 \, \mathrm{N} \, \mathrm{m}^{-1}$ |

material!(b; horizon=0.01, rho=2700, E=70e9, Gc=100)

# As loading condition for the specimen, a constant velocity of $0.6 \, \mathrm{m}\,\mathrm{s}^{-1}$ in $x$-direction is set for the bottom and top.
velocity_bc!(t -> -0.6, b, :bottom, 1)
velocity_bc!(t -> 0.6, b, :top, 1)

# Set the number of time steps for the Velocity Verlet algorithm to 500 time steps:
vv = VelocityVerlet(steps=500)

# Define the storage path:
root = joinpath(@__DIR__, "results", "tension_dynamic_fracture")

# Create the job:
job = Job(b, vv; path=root)

#-
# Submit the job to start simulations:
#md # ```julia
#md # submit(job)
#md # ```

# ![](../assets/tension_dynfrac.gif)
# (Visualization made with [ParaView](https://www.paraview.org))
