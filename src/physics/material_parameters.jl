@inline elasticity_parameters() = (:E, :nu, :G, :K, :λ, :μ)

@inline elasticity_kwargs() = (:E, :nu, :G, :K, :lambda, :mu)
@inline discretization_kwargs() = (:horizon, :rho)

"""
    material!(body, set_name; kwargs...)
    material!(body; kwargs...)

Assign material point parameters to points of `body`. If no `set_name` is specified, then
the parameters will be set for all points of the body.

# Arguments

- `body::AbstractBody`: [`Body`](@ref).
- `set_name::Symbol`: The name of a point set of this body.

# Keywords

Allowed keywords depend on the selected material model. Please look at the documentation
of the material you specified when creating the body.
The default material keywords are:

Material parameters:
- `horizon::Float64`: Radius of point interactions
- `rho::Float64`: Density
Elastic parameters:
- `E::Float64`: Young's modulus
- `nu::Float64`: Poisson's ratio
- `G::Float64`: Shear modulus
- `K::Float64`: Bulk modulus
- `lambda::Float64`: 1st Lamé parameter
- `mu::Float64`: 2nd Lamé parameter
Fracture parameters:
- `Gc::Float64`: Critical energy release rate
- `epsilon_c::Float64`: Critical strain

!!! note "Elastic parameters"
    Note that exactly two elastic parameters are required to specify a material.
    Please choose two out of the six allowed elastic parameters.

!!! note "Fracture parameters"
    To enable fracture in a simulation, define one of the allowed fracture parameters.
    If none are defined, fracture is disabled.

# Throws

- Error if a kwarg is not eligible for specification with the body material.

# Example

```julia-repl
julia> material!(body; horizon=3.0, E=2.1e5, rho=8e-6, Gc=2.7)

julia> body
1000-point Body{BBMaterial{NoCorrection}}:
  1 point set(s):
    1000-point set `all_points`
  1 point parameter(s):
    Parameters BBMaterial: δ=3.0, E=210000.0, nu=0.25, rho=8.0e-6, Gc=2.7
```
"""
function material! end

function material!(body::AbstractBody, set_name::Symbol; kwargs...)
    check_if_set_is_defined(body.point_sets, set_name)

    p = Dict{Symbol,Any}(kwargs)
    check_material_kwargs(body.mat, p)

    points = body.point_sets[set_name]
    params = get_point_params(body.mat, p)

    _material!(body, points, params)
    set_failure_permissions!(body, set_name, params)

    return nothing
end

function material!(body::AbstractBody; kwargs...)
    isempty(body.point_params) || empty!(body.point_params)

    material!(body, :all_points; kwargs...)

    return nothing
end

function _material!(b::AbstractBody, points::V, params::P) where {P,V}
    push!(b.point_params, params)
    id = length(b.point_params)
    b.params_map[points] .= id
    return nothing
end

function check_material_kwargs(mat::AbstractMaterial, p::Dict{Symbol,Any})
    allowed_kwargs = allowed_material_kwargs(mat)
    check_kwargs(p, allowed_kwargs)
    return nothing
end

function get_horizon(p::Dict{Symbol,Any})
    if !haskey(p, :horizon)
        throw(UndefKeywordError(:horizon))
    end
    δ::Float64 = float(p[:horizon])
    δ ≤ 0 && throw(ArgumentError("`horizon` should be larger than zero!\n"))
    return (; δ,)
end

function get_density(p::Dict{Symbol,Any})
    if !haskey(p, :rho)
        throw(UndefKeywordError(:rho))
    end
    rho::Float64 = float(p[:rho])
    rho ≤ 0 && throw(ArgumentError("`rho` should be larger than zero!\n"))
    return (; rho,)
end

function get_elastic_params(p::Dict{Symbol,Any})
    par = get_given_elastic_params(p)
    check_elastic_params(par)
    (; E, nu)= get_E_and_nu(par)

    G = E / (2 * (1 + nu))
    K = E / (3 * (1 - 2 * nu))
    λ = E * nu / ((1 + nu) * (1 - 2nu))
    μ = G
    # return named tuple containing all 6 elastic material parameters
    return (; E, nu, G, K, λ, μ)
end

function get_given_elastic_params(p::Dict{Symbol,Any})
    # get elastic parameters from dictionary
    if haskey(p, :E)
        E::Float64 = float(p[:E])
        E ≤ 0 && throw(ArgumentError("`E` should be larger than zero!\n"))
    else
        E = NaN
    end
    if haskey(p, :nu)
        nu::Float64 = float(p[:nu])
        nu ≤ 0 && throw(ArgumentError("`nu` should be larger than zero!\n"))
        nu ≥ 1 && throw(ArgumentError("too high value of `nu`! Condition: 0 < `nu` ≤ 1\n"))
    else
        nu = NaN
    end
    if haskey(p, :G)
        G::Float64 = float(p[:G])
        G ≤ 0 && throw(ArgumentError("`G` should be larger than zero!\n"))
    else
        G = NaN
    end
    if haskey(p, :K)
        K::Float64 = float(p[:K])
        K ≤ 0 && throw(ArgumentError("`K` should be larger than zero!\n"))
    else
        K = NaN
    end
    if haskey(p, :lambda)
        λ::Float64 = float(p[:lambda])
    else
        λ = NaN
    end
    if haskey(p, :mu)
        μ::Float64 = float(p[:mu])
        μ ≤ 0 && throw(ArgumentError("`μ` should be larger than zero!\n"))
    else
        μ = NaN
    end
    return (; E, nu, G, K, λ, μ)
end

function check_elastic_params(par)
    (; G, μ) = par
    # check if exactly 2 keywords out of {E, nu, G, K, λ, μ} are provided
        # Caution: μ & G are not independet parameters!
    if isfinite(G) && isfinite(μ)
        throw(ArgumentError("`G` and `μ` are defined! Please define either `G` or `μ`!"))
    elseif length(findall(isfinite, par)) < 2
        msg =  "Not enough material parameters defined!\n"
        msg *= "To characterize the material, two parameters are required!\n"
        throw(ArgumentError(msg))
    elseif length(findall(isfinite, par)) > 2
        msg =  "Too many material parameters defined!\n"
        msg *= "To characterize the material, only two parameters are required!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function get_E_and_nu(par)
    (; E, nu, G, K, λ, μ) = par
    # check which 2 parameters are provided & calculate E & nu

    if isfinite(E) && isfinite(nu)

    elseif isfinite(E) && isfinite(G)
        nu = E / (2 * G) - 1
    elseif isfinite(E) && isfinite(K)
        nu = (3 * K - E) / (6 * K)
    elseif isfinite(E) && isfinite(λ)
        nu = (-(E + λ) + sqrt((E + λ)^2 + 8 * λ^2)) / (4 * λ)
    elseif isfinite(E) && isfinite(μ)
        nu = E / (2 * μ) - 1
    elseif isfinite(nu) && isfinite(G)
        E = 2 * G * (1 + nu)
    elseif isfinite(nu) && isfinite(K)
        E = 3 * K * (1 - 2 * nu)
    elseif isfinite(nu) && isfinite(λ)
        E = (λ * (1 + nu) * (1 - 2 * nu)) / nu
    elseif isfinite(nu) && isfinite(μ)
        E = 2 * μ * (1 + nu)
    elseif isfinite(G) && isfinite(K)
        E = 9 * K * G / (3 * K + G)
        nu = (3 * K - 2 * G) / (2 * (3 * K + G))
    elseif isfinite(G) && isfinite(λ)
        E = G * (3 * λ + 2 * G) / (λ + G)
        nu = λ / (2 * (λ + G))
    elseif isfinite(K) && isfinite(λ)
        E = 9 * K * (K - λ) / (3 * K - λ)
        nu = λ / (3 * K - λ)
    elseif isfinite(K) && isfinite(μ)
        E = 9 * K * μ / (3 * K + μ)
        nu = (3 * K - 2 * μ) / (2 * (3 * K + μ))
    elseif isfinite(λ) && isfinite(μ)
        E = μ * (3 * λ + 2 * μ) / (λ + μ)
        nu = λ / (2 * (λ + μ))
    end

    # return named tuple containing E & nu
    return (; E, nu)
end

function log_material_parameters(param::AbstractPointParameters; indentation::Int=2)
    msg = msg_qty("horizon", param.δ; indentation=indentation)
    msg *= msg_qty("density", param.rho; indentation=indentation)
    msg *= msg_qty("Young's modulus", param.E; indentation=indentation)
    msg *= msg_qty("Poisson's ratio", param.nu; indentation=indentation)
    msg *= msg_qty("shear modulus", param.G; indentation=indentation)
    msg *= msg_qty("bulk modulus", param.K; indentation=indentation)
    return msg
end

function Base.show(io::IO, @nospecialize(params::AbstractPointParameters))
    print(io, "Parameters ", material_type(params), ": ")
    print(io, msg_fields_inline(params, (:δ, :E, :nu, :rho, :Gc)))
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain",
                   @nospecialize(params::AbstractPointParameters))
    if get(io, :compact, false)
        show(io, params)
    else
        println(io, typeof(params), ":")
        print(io, msg_fields(params))
    end
    return nothing
end
