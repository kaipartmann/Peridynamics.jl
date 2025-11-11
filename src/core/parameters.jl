function point_param_type(mat::AbstractMaterial)
    return throw(InterfaceError(mat, "point_param_type"))
end

function get_point_params(mat::AbstractMaterial, ::Dict{Symbol,Any})
    return throw(InterfaceError(mat, "get_point_params"))
end

function required_point_parameters(mat::Type{<:AbstractMaterial})
    return throw(InterfaceError(mat, "required_point_parameters"))
end

function allowed_material_kwargs(mat::AbstractMaterial)
    return throw(InterfaceError(mat, "allowed_material_kwargs"))
end

macro params(material, params)
    macrocheck_input_material(material)
    macrocheck_input_params(params)
    local _pre_checks = quote
        Peridynamics.typecheck_material($(esc(material)))
        Peridynamics.typecheck_params($(esc(material)), $(esc(params)))
    end
    local _pointparam_type = quote
        Peridynamics.point_param_type(::$(esc(material))) = $(esc(params))
    end
    local _get_pointparams = quote
        function Peridynamics.get_point_params(m::$(esc(material)), p::Dict{Symbol,Any})
            return $(esc(params))(m, p)
        end
    end
    local _post_checks = quote
        Peridynamics.constructor_check($(esc(material)), $(esc(params)))
    end
    exprs = Expr(:block, _pre_checks, _pointparam_type, _get_pointparams, _post_checks)
    return exprs
end

function macrocheck_input_material(material)
    material isa Symbol && return nothing
    (material isa Expr && material.head === :.) && return nothing
    (material isa Expr && material.head === :escape) && return nothing
    (material isa Expr && material.head === :curly) && return nothing
    return throw(ArgumentError("argument `$material` is not a valid material input!\n"))
end

function macrocheck_input_params(params)
    params isa Symbol && return nothing
    (params isa Expr && params.head === :.) && return nothing
    msg = "argument `$params` is not a valid point parameter input!\n"
    return throw(ArgumentError(msg))
end

function typecheck_material(::Type{Material}) where {Material}
    if !(Material <: AbstractMaterial)
        msg = "$Material is not a valid material type!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function typecheck_params(::Type{Material}, ::Type{Param}) where {Material,Param}
    if !(Param <: AbstractPointParameters)
        msg = "$Param is not a subtype of AbstractPointParameters!\n"
        throw(ArgumentError(msg))
    end
    parameters = fieldnames(Param)
    for req_param in required_point_parameters(Material)
        if !in(req_param, parameters)
            msg = "required parameter $req_param not found in $(Param)!\n"
            error(msg)
        end
    end
    return nothing
end

function constructor_check(::Type{Material}, ::Type{Param}) where {Material,Param}
    if !hasmethod(Param, Tuple{Material,Dict{Symbol,Any}})
        throw(InterfaceError(Param, "$(Param)(::$(Material), ::Dict{Symbol,Any})"))
    end
    return nothing
end

"""
    StandardPointParameters

$(internal_api_warning())

Type containing the material parameters for a standard peridynamics model using the
bond-based, ordinary state-based or non-ordinary state-based correspondence formulation of
peridynamics.

# Fields

- `δ::Float64`: Horizon.
- `rho::Float64`: Density.
- `E::Float64`: Young's modulus.
- `nu::Float64`: Poisson's ratio.
- `G::Float64`: Shear modulus.
- `K::Float64`: Bulk modulus.
- `λ::Float64`: 1st Lamé parameter.
- `μ::Float64`: 2nd Lamé parameter.
- `Gc::Float64`: Critical energy release rate.
- `εc::Float64`: Critical strain.
- `bc::Float64`: Bond constant.
"""
struct StandardPointParameters <: AbstractPointParameters
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

function StandardPointParameters(mat::AbstractMaterial, p::Dict{Symbol,Any})
    (; δ, rho, E, nu, G, K, λ, μ) = get_required_point_parameters(mat, p)
    (; Gc, εc) = get_frac_params(mat.dmgmodel, p, δ, K)
    # for a cubic neighborhood, the weighted volume is the integral
    # 8 δ^4 ∫_0^1 ∫_0^1 ∫_0^1 √(x² + y² + z²) dx dy dz
    # which results in 8 * δ^4 * 0.9605919564548167 (solved numerically)
    bc = 18 * K / (8 * δ^4 * 0.9605919564548167) # bond constant
    return StandardPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end
