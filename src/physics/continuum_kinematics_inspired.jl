"""
TODO
"""
struct CKIMaterial <: AbstractMaterial end

struct CKIPointParameters <: AbstractPointParameters
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
    C1::Float64
    C2::Float64
    C3::Float64
end

@inline point_param_type(::CKIMaterial) = CKIPointParameters

function get_point_params(::CKIMaterial, p::Dict{Symbol,Any})
    δ = get_horizon(p)
    rho = get_density(p)
    E, nu, G, K, λ, μ = get_elastic_params(p)
    Gc, εc = get_frac_params(p, δ, K)
    C1 = 30 / π * μ / δ^4
    C2 = 0.0
    C3 = 32 / π^4 * (λ - μ) / δ^12
    return CKIPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, C1, C2, C3)
end

@inline allowed_material_kwargs(::CKIMaterial) = DEFAULT_POINT_KWARGS
