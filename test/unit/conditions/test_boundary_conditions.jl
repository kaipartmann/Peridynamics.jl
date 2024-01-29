using Test
using Peridynamics

function create_bc(val::Float64)
    v1(t) = val * t
    Peridynamics.VelocityBC(v1, 1:10, 1)
end

bc1 = create_bc(1.0)
for t in [-1, 0, 1, Inf]
    @test bc1.fun(t) == 1.0 * t
end
bc2 = create_bc(2.0)
for t in [-1, 0, 1, Inf]
    @test bc1.fun(t) == 1.0 * t
    @test bc2.fun(t) == 2.0 * t
end
bc3 = Peridynamics.VelocityBC(1:10, 1) do t
    3.0 * t
end
bc4 = Peridynamics.VelocityBC(t -> 4.0t, 1:10, 1)
for t in [-1, 0, 1, Inf]
    @test bc1.fun(t) == 1.0 * t
    @test bc2.fun(t) == 2.0 * t
    @test bc3.fun(t) == 3.0 * t
    @test bc4.fun(t) == 4.0 * t
end

struct TestStorage <: Peridynamics.AbstractStorage
    velocity_half::Matrix{Float64}
    b_ext::Matrix{Float64}
end
s = TestStorage(zeros(3,10), zeros(3,10))
time = 0.1

Peridynamics.apply_bc!(s, bc, time)
@test s.velocity_half[1,1:10] == fill(v1(time), 10)
