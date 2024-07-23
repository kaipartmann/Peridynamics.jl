struct BANOSBMaterial <: AbstractBondAssociatedSystemMaterial end

struct BANOSBPointParameters <: AbstractPointParameters
    δ::Float64
    rho::Float64
    E::Float64
    nu::Float64
    G::Float64
    K::Float64
    λ::Float64
    μ::Float64
    Gc::Float64
    εc::Float64
    bc::Float64
end

function BANOSBPointParameters(::BANOSBMaterial, p::Dict{Symbol,Any})
    δ = get_horizon(p)
    rho = get_density(p)
    E, nu, G, K, λ, μ = get_elastic_params(p)
    Gc, εc = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return BANOSBPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params BANOSBMaterial BANOSBPointParameters
