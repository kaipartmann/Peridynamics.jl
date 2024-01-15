using Peridynamics, Test

dr1 = DynamicRelaxation(1)
@test dr1.n_steps == 1
@test dr1.Δt == 1
@test dr1.Λ == 1

dr2 = DynamicRelaxation(100, 1.5)
@test dr2.n_steps == 100
@test dr2.Δt == 1.5
@test dr2.Λ == 1

dr3 = DynamicRelaxation(100, 1.5; damping_factor=10)
@test dr3.n_steps == 100
@test dr3.Δt == 1.5
@test dr3.Λ == 10

dr4 = DynamicRelaxation(100; damping_factor=100)
@test dr4.n_steps == 100
@test dr4.Δt == 1
@test dr4.Λ == 100