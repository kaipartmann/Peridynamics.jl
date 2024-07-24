# # [Tensile test quasi-static](@id tutorial_tension_static)

#-
# Import the package:
using Peridynamics

#-
# Read and convert the [Abaqus FEM mesh of a tensile test](https://github.com/kaipartmann/Peridynamics.jl/blob/main/docs/src/assets/TensileTestMesh.inp)
# into a point cloud for the peridynamic model:
## insert your correct path to the downloaded mesh file!
inp_file = joinpath(@__DIR__, "..", "assets", "TensileTestMesh.inp");

# ![](https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/d2f333d4-ffa5-4227-aed9-542d49298c01)

#-
# Create a body with the points from the mesh: (The bond-based material model is used here.)
body = Body(BBMaterial(), inp_file)

# The element sets defined in Abaqus were converted into point sets of the Body:
point_sets(body)

# Specify the material parameters as:

# | material parameter | value |
# |:--------|:-------------|
# | Horizon $ δ $ | $0.01 \, \mathrm{m}$ |
# | Density $ρ$ | $ 2700 \,\mathrm{kg}\,\mathrm{m}^{-3}$ |
# | Young's modulus $E$ | $ 70 \cdot 10^{9} \, \mathrm{Pa}$ |
# | Griffith's parameter $G_c$ | $100 \, \mathrm{N} \, \mathrm{m}^{-1}$ |

material!(body; horizon=0.01, rho=2700, E=70e9, Gc=100)

#-
# As loading condition for the specimen, a constant force density of $1 \times 10^9 \, \mathrm{N}\,\mathrm{m}^{-3}$ in $x$-direction is set for the bottom and top.
forcedensity_bc!(t -> -3e11, body, :bottom, 1)
forcedensity_bc!(t -> 3e11, body, :top, 1)

#-
# Do not allow failure in the entire body:
failure_permit!(body, false)

#-
# We set the number of time steps for the dynamic relaxation algorithm to 500 time steps.
dr = DynamicRelaxation(steps=500, damping_factor=0.2)

# Create the job:
job = Job(body, dr; path="results/tension_static")

#-
#md # ```julia
#md # submit(job)
#md # ```

# ![](https://github.com/kaipartmann/Peridynamics.jl/assets/68582683/d74d60fb-c792-4d5b-a683-9448d6c92e83)
