@inline elasticity_parameters() = (:E, :nu, :G, :K, :λ, :μ)

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
- `horizon::Float64`: Radius of point interactions.
- `rho::Float64`: Density.
Elastic parameters:
- `E::Float64`: Young's modulus.
- `nu::Float64`: Poisson's ratio.
- `G::Float64`: Shear modulus.
- `K::Float64`: Bulk modulus.
- `lambda::Float64`: 1st Lamé parameter.
- `mu::Float64`: 2nd Lamé parameter.
Fracture parameters:
- `Gc::Float64`: Critical energy release rate.
- `epsilon_c::Float64`: Critical strain.

!!! note "Elastic parameters"
    Note that exactly two elastic parameters are required to specify a material.

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
    check_model_params(body.mat)
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
    check_kwargs(p, all_material_kwargs(mat))
    return nothing
end

"""
    all_material_kwargs(mat)

$(internal_api_warning())

Every keyword `material!` accepts for a material: the keywords its own parameter
declarations consume, plus the keywords of its constitutive model and of its damage model.
A keyword is accepted exactly when something reads it, and a keyword two of them declare is
an error naming both owners.
"""
function all_material_kwargs(mat::AbstractMaterial)
    material = allowed_material_kwargs(mat)
    cm = constitutive_param_kwargs(get_constitutive_model(mat))
    dmg = damage_param_kwargs(get_dmgmodel(mat))
    check_kwarg_owners(mat, material, cm, dmg)
    return (material..., cm..., dmg...)
end

function check_kwarg_owners(mat, material, cm, dmg)
    owners = Dict{Symbol,Vector{String}}()
    for (kwargs, owner) in ((material, "the material `$(nameof(typeof(mat)))`"),
                            (cm, "the constitutive model"),
                            (dmg, "the damage model"))
        for kwarg in kwargs
            push!(get!(Vector{String}, owners, kwarg), owner)
        end
    end
    for (kwarg, names) in owners
        length(names) > 1 || continue
        msg = "the `material!` keyword `$(kwarg)` is declared more than once: by "
        msg *= "$(join(names, " and "))!\n"
        msg *= "  Rename the keyword with `@kwarg` in one of the declarations.\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function get_horizon(; horizon)
    δ::Float64 = float(horizon)
    δ ≤ 0 && throw(ArgumentError("`horizon` should be larger than zero!\n"))
    return (; δ,)
end

function get_density(; rho)
    ρ::Float64 = float(rho)
    ρ ≤ 0 && throw(ArgumentError("`rho` should be larger than zero!\n"))
    return (; rho=ρ)
end

"""
    DiscretizationParameters

$(extension_api_note())

Parameter block of the horizon `δ` and the density `rho`, the two parameters every
peridynamic material needs. See [`@params_fields`](@ref).

$(block_table(DiscretizationParameters))
"""
@params_fields DiscretizationParameters begin
    @derived (; δ, rho) = get_discretization_params(; horizon, rho)
    @log "horizon" δ
    @log "density" rho
end

function get_discretization_params(; horizon, rho)
    return (; get_horizon(; horizon)..., get_density(; rho)...)
end

"""
    ElasticParameters

$(extension_api_note())

Parameter block of the six elastic parameters. They are resolved together, because any two
of the six keywords `E`, `nu`, `G`, `K`, `lambda` and `mu` determine all of them. See
[`@params_fields`](@ref).

$(block_table(ElasticParameters))
"""
@params_fields ElasticParameters begin
    @derived (; E, nu, G, K, λ, μ) = get_elastic_params(; E, nu, G, K, lambda, mu)
    @log "Young's modulus" E
    @log "Poisson's ratio" nu
    @log "shear modulus" G
    @log "bulk modulus" K
end

function get_elastic_params(; E=nothing, nu=nothing, G=nothing, K=nothing, lambda=nothing,
                            mu=nothing)
    return resolve_elastic_params(get_given_elastic_params(; E, nu, G, K, lambda, mu))
end

# the six elastic parameters that follow from the two that were given
function resolve_elastic_params(given)
    check_elastic_params(given)
    (; E, nu) = get_E_and_nu(given)

    G = E / (2 * (1 + nu))
    K = E / (3 * (1 - 2 * nu))
    λ = E * nu / ((1 + nu) * (1 - 2nu))
    μ = G
    # return named tuple containing all 6 elastic material parameters
    return (; E, nu, G, K, λ, μ)
end

function get_given_elastic_params(; E=nothing, nu=nothing, G=nothing, K=nothing,
                                  lambda=nothing, mu=nothing)
    # a keyword that was not specified is `nothing` and becomes `NaN`, so that the number of
    # given parameters is `count(isfinite, ...)`
    _E = given_elastic_param(E, "E")
    _nu = given_elastic_param(nu, "nu")
    _nu ≥ 1 && throw(ArgumentError("too high value of `nu`! Condition: 0 < `nu` ≤ 1\n"))
    _G = given_elastic_param(G, "G")
    _K = given_elastic_param(K, "K")
    # the 1st Lamé parameter is the only one that may be negative
    λ::Float64 = isnothing(lambda) ? NaN : float(lambda)
    _μ = given_elastic_param(mu, "μ")
    return (; E=_E, nu=_nu, G=_G, K=_K, λ, μ=_μ)
end

function given_elastic_param(value, name::String)
    isnothing(value) && return NaN
    x::Float64 = float(value)
    x ≤ 0 && throw(ArgumentError("`$(name)` should be larger than zero!\n"))
    return x
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

function Base.show(io::IO, @nospecialize(params::AbstractPointParameters))
    print(io, nameof(typeof(params)), ": ")
    props = Tuple(s for s in (:δ, :E, :nu, :rho, :Gc) if hasproperty(params, s))
    values = Tuple(getproperty(params, s) for s in props)
    print(io, _msg_fields_inline(props, values))
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain",
                   @nospecialize(params::AbstractPointParameters))
    if get(io, :compact, false)
        show(io, params)
    else
        println(io, nameof(typeof(params)), ":")
        props = flat_param_property_names(params)
        values = Tuple(getproperty(params, s) for s in props)
        print(io, _msg_fields(props, values))
    end
    return nothing
end
