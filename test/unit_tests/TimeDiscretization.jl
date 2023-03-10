using Peridynamics, Test

td1 = TimeDiscretization(1)
@test td1.n_timesteps == 1
@test td1.Δt == -1
@test td1.alg == :verlet

td2 = TimeDiscretization(2; alg=:dynrelax)
@test td2.n_timesteps == 2
@test td2.Δt == 1
@test td2.alg == :dynrelax

@test_throws DomainError TimeDiscretization(2; alg = :somethingwrong)

td3 = TimeDiscretization(3, 0.1)
@test td3.n_timesteps == 3
@test td3.Δt == 0.1
@test td3.alg == :verlet

msg = "for dynamic relaxation a time step of Δt = 1 is recommended!"
td4 = @test_warn msg TimeDiscretization(3, 0.2; alg=:dynrelax)
@test td4.n_timesteps == 3
@test td4.Δt == 0.2
@test td4.alg == :dynrelax

@test_throws DomainError TimeDiscretization(4, 0.4; alg = :somethingwrong)
