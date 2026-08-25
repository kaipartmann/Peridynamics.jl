@testitem "get_horizon" begin
    (; δ) = Peridynamics.get_horizon(; horizon=1)
    @test δ ≈ 1.0

    # a keyword that was not given is not forwarded, so the provider reports it
    @test_throws UndefKeywordError(:horizon) Peridynamics.get_horizon()

    msg = "`horizon` should be larger than zero!\n"
    @test_throws ArgumentError(msg) Peridynamics.get_horizon(; horizon=0)
end

@testitem "get_density" begin
    (; rho) = Peridynamics.get_density(; rho=1)
    @test rho ≈ 1.0

    @test_throws UndefKeywordError(:rho) Peridynamics.get_density()

    msg = "`rho` should be larger than zero!\n"
    @test_throws ArgumentError(msg) Peridynamics.get_density(; rho=0)
end

@testitem "get_elastic_params" begin
    p = (; E=1, nu=0.25)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; E=1, G=0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; E=1, K=2/3)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; E=1, lambda=0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; E=1, mu=0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; nu=0.25, G=0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; nu=0.25, K=2/3)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; nu=0.25, lambda=0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; nu=0.25, mu=0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; G=0.4, K=2/3)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; G=0.4, lambda=0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; G=0.4, mu=0.4)
    @test_throws ArgumentError Peridynamics.get_elastic_params(; p...)

    p = (; K=2/3, lambda=0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; K=2/3, mu=0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; lambda=0.4, mu=0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; E=1, nu=0.25, G=NaN, K=NaN, lambda=NaN, mu=NaN)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = (; E=1, nu=0.25, G=NaN, K=3, lambda=NaN, mu=NaN)
    @test_throws ArgumentError Peridynamics.get_elastic_params(; p...)

    p = (; E=1, nu=NaN, G=NaN, K=NaN, lambda=NaN, mu=NaN)
    @test_throws ArgumentError Peridynamics.get_elastic_params(; p...)

    p = (; E=210e9, nu=0.25)
    E, nu, G, K, λ, μ = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 2.1e11
    @test nu ≈ 0.25
    @test G ≈ 8.4e10
    @test K ≈ 1.4e11
    @test λ ≈ 8.4e10
    @test μ ≈ 8.4e10

    p = (; E=27e9, nu=0.2)
    E, nu, G, K, λ, μ = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 2.7e10
    @test nu ≈ 0.2
    @test G ≈ 1.125e10
    @test K ≈ 1.5e10
    @test λ ≈ 7.5e9
    @test μ ≈ 1.125e10

    p = (; G=1, K=1)
    E, nu, G, K, λ, μ = Peridynamics.get_elastic_params(; p...)
    @test E ≈ 2.25
    @test nu ≈ 0.125
    @test G ≈ 1.0
    @test K ≈ 1.0
    @test λ ≈ 1/3
    @test μ ≈ 1.0

    p = (; E=0, nu=0.25)
    msg = "`E` should be larger than zero!\n"
    @test_throws ArgumentError Peridynamics.get_elastic_params(; p...)

    p = (; E=1, nu=0)
    msg = "`nu` should be larger than zero!\n"
    @test_throws ArgumentError Peridynamics.get_elastic_params(; p...)

    p = (; E=1, nu=1.1)
    msg = "too high value of `nu`! Condition: 0 < `nu` ≤ 1\n"
    @test_throws ArgumentError Peridynamics.get_elastic_params(; p...)

    p = (; E=1, G=0)
    msg = "`G` should be larger than zero!\n"
    @test_throws ArgumentError Peridynamics.get_elastic_params(; p...)

    p = (; E=1, K=0)
    msg = "`K` should be larger than zero!\n"
    @test_throws ArgumentError Peridynamics.get_elastic_params(; p...)

    p = (; E=1, mu=0)
    msg = "`μ` should be larger than zero!\n"
    @test_throws ArgumentError Peridynamics.get_elastic_params(; p...)
end

@testitem "log_param_property" begin
    p = Dict{Symbol,Any}(:E => 1, :nu => 0.25, :rho => 1, :horizon => 1)
    param = Peridynamics.get_point_params(BBMaterial(), p)
    msg = Peridynamics.log_param_property(Val(:randomthing), param; indentation=0)
    @test msg == ""
    msg = Peridynamics.log_param_property(Val(:δ), param; indentation=0)
    @test contains(msg, "horizon")
    @test contains(msg, "1")
    msg = Peridynamics.log_param_property(Val(:E), param; indentation=0)
    @test contains(msg, "Young's modulus")
    @test contains(msg, "1")
end

@testitem "show point parameters" setup=[Fixtures] begin
    body = Fixtures.tetra4(OSBMaterial(); horizon=2.0, rho=1.0, E=1.0, nu=0.25)
    params = body.point_params[1]
    name = string(nameof(typeof(params)))
    msg = sprint(show, params)
    @test contains(msg, name * ": ") && !contains(msg, "\n")
    @test contains(msg, "δ=2.0") && contains(msg, "E=1.0") && contains(msg, "nu=0.25")
    msg = sprint(show, MIME("text/plain"), params)
    @test contains(msg, name * ":\n")
    @test contains(msg, "δ") && contains(msg, "rho") && contains(msg, "bc")
    msg = sprint(show, MIME("text/plain"), params; context=:compact => true)
    @test contains(msg, name * ": ") && !contains(msg, "\n")
end
