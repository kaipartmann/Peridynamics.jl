# # [Brazilian Test](@id tutorial_brazilian_test)

# This tutorial sets up the Brazilian Test experiment, commonly used to investigate fracture
# of brittle materials like ultra-high performance concrete.
# Therefore a cylindrical specimen is
# loaded by two opposing forces applying pressure on the cross section of the specimen.

# To start, we import the package.
using Peridynamics

# Then we define the geometrical parameters of the specimen which are the diameter `Ø` and
# the length `l` as well as the point spacing `Δx` of the model.
# The parameter `b` is used later defining point sets for the boundary conditions.
Ø = 0.05 # [m]
l = 0.015 # [m]
Δx = Ø/61 # [m]
b = 0.017 # [m]

# With the diameter, length and point spacing, we create the cylindrical body.
pos, vol = uniform_cylinder(Ø, l, Δx)
cyl = Body(BBMaterial(), pos, vol)

# The horizon is specified in relation to the point spacing.
δ = 3.015Δx
# Then the material parameters are set.
material!(cyl; horizon=δ, E=50e9, rho=2400, Gc=140)

# To apply the opposing forces, two point sets are generated.
point_set!(p -> p[1] ≥ Ø/2-3Δx && abs(p[2]) ≤ 1.1*b/2, cyl, :set_1)
point_set!(p -> p[1] ≤ -Ø/2+3Δx && abs(p[2]) ≤ 1.1*b/2, cyl, :set_2)

# The functions for the velocity boundary conditions, that have the same value but act in
# different directions, are defined next.
v_set1(p, t) = - 2 * exp(-t/0.00002) * (-1/0.02^2 * p[2]^2 + 1)
v_set2(p, t) =  2 * exp(-t/0.00002) * (-1/0.02^2 * p[2]^2 + 1)

# Then these functions are set as the boundary conditions in x-direction, while the
# velocites of the points in y- and z-direction are 0.
velocity_bc!(v_set1, cyl, :set_1, :x)
velocity_bc!(t -> 0, cyl, :set_1, :y)
velocity_bc!(t -> 0, cyl, :set_1, :z)
velocity_bc!(v_set2, cyl, :set_2, :x)
velocity_bc!(t -> 0, cyl, :set_2, :y)
velocity_bc!(t -> 0, cyl, :set_2, :z)

# The Velocity Verlet algotihm is employed as time integration method where 6000 time steps
# are calculated.
vv = VelocityVerlet(steps=6000)

# Finally the job is defined and submitted.

#md # ```julia
#md # job = Job(cyl, vv; path="results/brazilian_bb_uniform")
#md # submit(job)
#md # ```
