# # [Fragmenting Cylinder](@id tutorial_cylinder)

# A cylinder fragmenting into many pieces, inspired by the
#  [peridigm](https://github.com/peridigm/peridigm/tree/master) example.

# We use the point cloud created by [peridigm](https://github.com/peridigm/peridigm/tree/master)
# for their corresponding example.
# You can download the file containing the data from
# [this page](https://github.com/peridigm/peridigm/blob/master/examples/fragmenting_cylinder/fragmenting_cylinder.txt)
# on their repository.

# To start, import the package:
using Peridynamics

# First a function is written that can read our `.txt`-file containing the data of the
# point cloud and convert it into the position and volume, which we need to define our
# `Body`.
# For this we need the package `DelimitedFiles.jl`.
using DelimitedFiles
function fragmenting_cylinder_geometry(input_mesh_file::AbstractString)
    input_raw = readdlm(input_mesh_file)
    position = copy(input_raw[:, 1:3]')
    volume = copy(input_raw[:, 5])
    return position, volume
end

# Now we specify the storage path of the file.
input_mesh_file = joinpath(@__DIR__, "..", "assets", "fragmenting_cylinder.txt")
# To get the information we need, we use our function defined above.
position, volume = fragmenting_cylinder_geometry(input_mesh_file)
# Using this data, we can create a `Body` which represents the cylinder.
body = Body(BBMaterial(), position, volume)
# We specify the material parameters of the cylinder.
material!(body, horizon=0.00417462, rho=7800, E=195e9, epsilon_c=0.02)
# Then some initial velocity conditions in x-, y- and z-direction are employed to provoke
# the fracture of the cylinder.
velocity_ic!(p -> (200-50*((p[3]/0.05)-1)^2)*cos(atan(p[2],p[1])), body, :all_points, :x)
velocity_ic!(p -> (200-50*((p[3]/0.05)-1)^2)*sin(atan(p[2],p[1])), body, :all_points, :y)
velocity_ic!(p -> 100*((p[3]/0.05)-1), body, :all_points, :z)
# We employ the Velocity Verlet algorithm for a total time span of 2.5e-4 seconds.
vv = VelocityVerlet(time=2.5e-4)
# Finally the job is created
job = Job(body, vv; path="results/fragmenting_cylinder", freq=10)
# and subsequently submitted.
#md # ```julia
#md # submit(job)
#md # ```
