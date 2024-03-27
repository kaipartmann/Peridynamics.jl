const ELASTIC_KWARGS = (:E, :nu)
const FRAC_KWARGS = (:Gc, :epsilon_c)
const DEFAULT_POINT_KWARGS = (:horizon, :rho, ELASTIC_KWARGS..., FRAC_KWARGS...)

function get_horizon(p::Dict{Symbol,Any})
    if !haskey(p, :horizon)
        throw(UndefKeywordError(:horizon))
    end
    δ::Float64 = float(p[:horizon])
    δ ≤ 0 && throw(ArgumentError("`horizon` should be larger than zero!\n"))
    return δ
end

function get_density(p::Dict{Symbol,Any})
    if !haskey(p, :rho)
        throw(UndefKeywordError(:rho))
    end
    rho::Float64 = float(p[:rho])
    rho ≤ 0 && throw(ArgumentError("`rho` should be larger than zero!\n"))
    return rho
end

function get_elastic_params(p::Dict{Symbol,Any})
    if !haskey(p, :E) || !haskey(p, :nu)
        msg = "insufficient keywords for calculation of elastic parameters!\n"
        msg *= "The keywords `E` (elastic modulus) and `nu` (poisson ratio) are needed!\n"
        throw(ArgumentError(msg))
    end
    E::Float64 = float(p[:E])
    E ≤ 0 && throw(ArgumentError("`E` should be larger than zero!\n"))
    nu::Float64 = float(p[:nu])
    nu ≤ 0 && throw(ArgumentError("`nu` should be larger than zero!\n"))
    nu ≥ 1 && throw(ArgumentError("too high value of `nu`! Condition: 0 < `nu` ≤ 1\n"))
    G = E / (2 * (1 + nu))
    K = E / (3 * (1 - 2 * nu))
    λ = E * nu / ((1 + nu) * (1 - 2nu))
    μ = G
    return E, nu, G, K, λ, μ
end
