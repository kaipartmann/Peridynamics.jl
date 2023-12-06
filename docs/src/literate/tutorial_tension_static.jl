# # [Mode I tension quasi-static](@id tutorial_tension_static)

#-
# Import the package:
using Peridynamics

#-
# Read and convert the [Abaqus FEM mesh of a tensile test](https://github.com/kaipartmann/Peridynamics.jl/blob/main/docs/src/assets/TensileTestMesh.inp):
## insert your correct path to the downloaded mesh file!
pc = read_inp(joinpath(@__DIR__, "..", "assets", "TensileTestMesh.inp"))
# ![](../assets/TensileTestMesh.png)

#-
# Do not allow failure in the entire point cloud:
pc.failure_flag .= false

#-
# Define a bond-based material with
# - Horizon $\delta = 0.01\,\mathrm{m}$
# - Density $\rho = 2700\,\mathrm{kg}\,\mathrm{m}^{-3}$
# - Youngs modulus $E = 70 \times 10^9 \, \mathrm{Pa}$
# - Griffith's parameter $G_c = 100 \, \mathrm{N} \, \mathrm{m}^{-1}$
mat = BBMaterial(; horizon=0.01, rho=2700, E=70e9, Gc=100)

#-
# As loading condition for the specimen, a constant force density of $1 \times 10^9 \, \mathrm{N}\,\mathrm{m}^{-3}$ in $x$-direction is set for the bottom and top.
# Note, that element sets defined in Abaqus are converted to `point_sets` of the `PointCloud`.
bcs = [
    ForceDensityBC(t -> -1e9, pc.point_sets["bottom"], 1),
    ForceDensityBC(t -> 1e9, pc.point_sets["top"], 1),
]

#-
# We set the number of time steps for the dynamic relaxation algorithm to 500 time steps.
dr = DynamicRelaxation(500)

#-
# The results of our analysis should be saved in the directory `"results/TensileTestStatic"` every 10'th time step.
name = "TensileTestStatic"
path = joinpath("results", name)
!ispath(path) && mkpath(path) # create the path if it does not exist
es = ExportSettings(path, 10)

#-
# Run a single body analysis:
#-
job = PDSingleBodyAnalysis(name=name, pc=pc, mat=mat, bcs=bcs, td=dr, es=es)

#-
#md # ```julia
#md # submit(job)
#md # ```

# ![](../assets/tension_static.gif)
# (Visualization made with [ParaView](https://www.paraview.org))
