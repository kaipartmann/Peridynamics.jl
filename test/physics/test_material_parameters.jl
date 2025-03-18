@testitem "get_horizon" begin
    p = Dict{Symbol,Any}(:horizon => 1)
    (; δ) = Peridynamics.get_horizon(p)
    @test δ ≈ 1.0

    p = Dict{Symbol,Any}(:nothorizon => 1)
    @test_throws UndefKeywordError(:horizon) Peridynamics.get_horizon(p)

    p = Dict{Symbol,Any}(:horizon => 0)
    msg = "`horizon` should be larger than zero!\n"
    @test_throws ArgumentError(msg) Peridynamics.get_horizon(p)
end

@testitem "get_density" begin
    p = Dict{Symbol,Any}(:rho => 1)
    (; rho) = Peridynamics.get_density(p)
    @test rho ≈ 1.0

    p = Dict{Symbol,Any}(:notrho => 1)
    @test_throws UndefKeywordError(:rho) Peridynamics.get_density(p)

    p = Dict{Symbol,Any}(:rho => 0)
    msg = "`rho` should be larger than zero!\n"
    @test_throws ArgumentError(msg) Peridynamics.get_density(p)
end

@testitem "get_elastic_params" begin
    p = Dict{Symbol,Any}(:E => 1, :nu => 0.25)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:E => 1, :G => 0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:E => 1, :K => 2/3)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:E => 1, :lambda => 0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:E => 1, :mu => 0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:nu => 0.25, :G => 0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:nu => 0.25, :K => 2/3)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:nu => 0.25, :lambda => 0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:nu => 0.25, :mu => 0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:G => 0.4, :K => 2/3)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:G => 0.4, :lambda => 0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:G => 0.4, :mu => 0.4)
    @test_throws ArgumentError Peridynamics.get_elastic_params(p)

    p = Dict{Symbol,Any}(:K => 2/3, :lambda => 0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:K => 2/3, :mu => 0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:lambda => 0.4, :mu => 0.4)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:E => 1, :nu => 0.25, :G => NaN, :K => NaN, :lambda => NaN, :mu => NaN)
    (; E, nu, G, K, λ, μ) = Peridynamics.get_elastic_params(p)
    @test E ≈ 1.0
    @test nu ≈ 0.25
    @test G ≈ 0.4
    @test K ≈ 2/3
    @test λ ≈ 0.4
    @test μ ≈ 0.4

    p = Dict{Symbol,Any}(:E => 1, :nu => 0.25, :G => NaN, :K => 3, :lambda => NaN, :mu => NaN)
    @test_throws ArgumentError Peridynamics.get_elastic_params(p)

    p = Dict{Symbol,Any}(:E => 1, :nu => NaN, :G => NaN, :K => NaN, :lambda => NaN, :mu => NaN)
    @test_throws ArgumentError Peridynamics.get_elastic_params(p)

    p = Dict{Symbol,Any}(:E => 210e9, :nu => 0.25)
    E, nu, G, K, λ, μ = Peridynamics.get_elastic_params(p)
    @test E ≈ 2.1e11
    @test nu ≈ 0.25
    @test G ≈ 8.4e10
    @test K ≈ 1.4e11
    @test λ ≈ 8.4e10
    @test μ ≈ 8.4e10

    p = Dict{Symbol,Any}(:E => 27e9, :nu => 0.2)
    E, nu, G, K, λ, μ = Peridynamics.get_elastic_params(p)
    @test E ≈ 2.7e10
    @test nu ≈ 0.2
    @test G ≈ 1.125e10
    @test K ≈ 1.5e10
    @test λ ≈ 7.5e9
    @test μ ≈ 1.125e10

    p = Dict{Symbol,Any}(:G => 1, :K => 1)
    E, nu, G, K, λ, μ = Peridynamics.get_elastic_params(p)
    @test E ≈ 2.25
    @test nu ≈ 0.125
    @test G ≈ 1.0
    @test K ≈ 1.0
    @test λ ≈ 1/3
    @test μ ≈ 1.0

    p = Dict{Symbol,Any}(:E => 0, :nu => 0.25)
    msg = "`E` should be larger than zero!\n"
    @test_throws ArgumentError Peridynamics.get_elastic_params(p)

    p = Dict{Symbol,Any}(:E => 1, :nu => 0)
    msg = "`nu` should be larger than zero!\n"
    @test_throws ArgumentError Peridynamics.get_elastic_params(p)

    p = Dict{Symbol,Any}(:E => 1, :nu => 1.1)
    msg = "too high value of `nu`! Condition: 0 < `nu` ≤ 1\n"
    @test_throws ArgumentError Peridynamics.get_elastic_params(p)

    p = Dict{Symbol,Any}(:E => 1, :G => 0)
    msg = "`G` should be larger than zero!\n"
    @test_throws ArgumentError Peridynamics.get_elastic_params(p)

    p = Dict{Symbol,Any}(:E => 1, :K => 0)
    msg = "`K` should be larger than zero!\n"
    @test_throws ArgumentError Peridynamics.get_elastic_params(p)

    p = Dict{Symbol,Any}(:E => 1, :mu => 0)
    msg = "`μ` should be larger than zero!\n"
    @test_throws ArgumentError Peridynamics.get_elastic_params(p)
end
