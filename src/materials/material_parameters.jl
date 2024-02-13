function get_horizon(p::Dict{Symbol,Any})
    if !haskey(p, :horizon)
        throw(UndefKeywordError(:horizon))
    end
    δ::Float64 = float(p[:horizon])
    return δ
end

function get_density(p::Dict{Symbol,Any})
    if !haskey(p, :rho)
        throw(UndefKeywordError(:rho))
    end
    rho::Float64 = float(p[:rho])
    return rho
end

function get_elastic_params(p::Dict{Symbol,Any})
    local E::Float64
    local nu::Float64
    local G::Float64
    local K::Float64
    local λ::Float64
    local μ::Float64

    if haskey(p, :E) && haskey(p, :nu)
        E = float(p[:E])
        nu = float(p[:nu])
        G = E / (2 * (1 + nu))
        K = E / (3 * (1 - 2 * nu))
        λ = E * nu / ((1 + nu) * (1 - 2nu))
        μ = G
    elseif haskey(p, :G) && haskey(p, :K)
        msg = "parameter calculation with keywords G and K not yet implemented!\n"
        msg *= "The following combinations of keywords are possible:\n"
        msg *= "  * `E` (elastic modulus) & `nu` (poisson ratio)\n"
        throw(ArgumentError(msg))
    elseif haskey(p, :lambda) && haskey(p, :mu)
        msg = "parameter calculation with keywords lambda and mu not yet implemented!\n"
        msg *= "The following combinations of keywords are possible:\n"
        msg *= "  * `E` (elastic modulus) & `nu` (poisson ratio)\n"
        throw(ArgumentError(msg))
    else
        msg = "insufficient keywords for calculation of elastic parameters!\n"
        msg *= "The following combinations of keywords are possible:\n"
        msg *= "  * `E` (elastic modulus) & `nu` (poisson ratio)\n"
        # msg *= "  * `G` (shear modulus) & `K` (bulk modulus)\n"
        # msg *= "  * `lambda` (1. Lamé parameter) & `mu` (2. Lamé parameter)\n"
        throw(ArgumentError(msg))
    end

    return E, nu, G, K, λ, μ
end
