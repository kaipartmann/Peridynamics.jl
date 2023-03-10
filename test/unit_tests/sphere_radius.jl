using Peridynamics, Test

# volume of sphere with radius r
V(r) = 4/3 * π * r^3

# testing values
r_vals = (0, 1, π, 100 * rand())
for r in r_vals
    @test Peridynamics.sphere_radius(V(r)) ≈ r
end

# error for negative values
# testing values
r_vals = (-1, -π, -100 * rand())
for r in r_vals
    @test_throws DomainError Peridynamics.sphere_radius(V(r))
end
