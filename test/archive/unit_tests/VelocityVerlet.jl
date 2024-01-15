using Peridynamics, Test

vv1 = VelocityVerlet(1)
@test vv1.n_steps == 1
@test vv1.Δt == -1
@test vv1.safety_factor == 0.7

vv2 = VelocityVerlet(100, 1.5)
@test vv2.n_steps == 100
@test vv2.Δt == 1.5
@test vv2.safety_factor == 0.7

vv3 = VelocityVerlet(100, 1.5; safety_factor=0.9)
@test vv3.n_steps == 100
@test vv3.Δt == 1.5
@test vv3.safety_factor == 0.9

vv4 = VelocityVerlet(100; safety_factor=0.9)
@test vv4.n_steps == 100
@test vv4.Δt == -1
@test vv4.safety_factor == 0.9