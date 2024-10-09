function point_param_type(mat::AbstractMaterial)
    throw(MethodError(point_param_type, mat))
end

function material_type(params::AbstractPointParameters)
    throw(MethodError(material_type, params))
end

function get_point_params(mat::AbstractMaterial, ::Dict{Symbol,Any})
    throw(MethodError(get_point_params, mat))
end

function extract_parameters(mat::AbstractMaterial,
                            custom_params::Vector{Tuple{Symbol,DataType}},
                            p::Dict{Symbol,Any})
    required_params = extract_required_parameters(mat, p)
    params = (; required_params...)
    for (name, type) in custom_params
        value = parse_parameter(mat, name, type, p)
        if isnothing(value)
            value = calc_parameter(mat, Val(name), required_params)
        end
        check_parameter(mat, Val(name), value)
        params = (; params..., name => value)
    end
    return params
end

function set_kw(::Type, ::Val{keyword}) where {keyword}
    return keyword
end

function extract_required_parameters(::AbstractMaterial, p::Dict{Symbol,Any})
    params = (; get_horizon(p)..., get_density(p)..., get_elastic_params(p)...)
    params = (; params..., get_frac_params(p, params)...)
    return params
end

function parse_parameter(::M, name::Symbol, ::Type{T}, p::Dict{Symbol,Any}) where {M,T}
    keyword = set_kw(M, Val(name))
    isnothing(keyword) && return nothing
    haskey(p, keyword) || return nothing
    value::T = convert(T, p[keyword])
    return value
end

function calc_parameter(::AbstractMaterial, ::Val{P}, elparams) where {P}
    msg = "point parameter `$P` is not specified and has no calculation method!\n"
    return throw(ArgumentError(msg))
end

function check_parameter(::AbstractMaterial, ::Val{P}, value) where {P}
    return nothing
end

function required_point_parameters(::Val{M}) where {M}
    required_params = (:(δ::Float64), # horizon
                       :(rho::Float64), # density
                       :(E::Float64), # Young's modulus
                       :(nu::Float64), # Poisson's ratio
                       :(G::Float64), # shear modulus
                       :(K::Float64), # bulk modulus
                       :(λ::Float64), # Lamé's first parameter
                       :(μ::Float64), # Lamé's second parameter
                       :(Gc::Float64), # fracture energy
                       :(εc::Float64))
    return required_params
end

function required_point_parameter_names(mat::Val{M}) where {M}
    required_params = required_point_parameters(mat)
    param_names = Tuple(extract_parameter_name(field) for field in required_params)
    return param_names
end

function required_point_parameter_types(mat::Val{M}) where {M}
    required_params = required_point_parameters(mat)
    param_types = Tuple(extract_parameter_type(field) for field in required_params)
    return param_types
end

function extract_parameter_name(expr::Expr)::Symbol
    if expr.head !== :(::)
        throw(ArgumentError("invalid expression for parameter extraction!\n"))
    end
    return expr.args[1]
end

function extract_parameter_type(expr::Expr)::DataType
    if expr.head !== :(::)
        throw(ArgumentError("invalid expression for parameter extraction!\n"))
    end
    return eval(expr.args[2])
end

function expand_point_parameters(struct_expr, matval)
    # get the struct name and the input fields
    struct_header = struct_expr.args[2]
    input_fields = filter(field -> field isa Expr, struct_expr.args[3].args)
    # construct the fields block
    required_params = required_point_parameters(matval)
    required_param_names = required_point_parameter_names(matval)
    fields_expr = Vector{Expr}()
    custom_fields = Vector{Tuple{Symbol,DataType}}()
    for field in input_fields
        if field.head === :(::)
            push!(fields_expr, field)
            name = extract_parameter_name(field)
            if !in(name, required_param_names)
                custom_field = (name, extract_parameter_type(field))
                push!(custom_fields, custom_field)
            end
        elseif field.head === :macrocall && field.args[1] === Symbol("@required_params")
            push!(fields_expr, required_params...)
        else
            msg = "unexpected or untyped field: $field\n"
            throw(ArgumentError(msg))
        end
    end
    fields_block = Expr(:block, fields_expr...)
    # construct the expanded struct
    _struct = Expr(:struct, struct_expr.args[1], struct_header, fields_block)
    resulting_struct = :(Base.@kwdef $(_struct))
    return resulting_struct, custom_fields
end

macro required_params()
    return nothing
end

macro params(material, params_expr)
    macrocheck_input_material(material)
    macrocheck_input_params(params_expr)
    params = extract_struct_name(params_expr)
    matval = extract_material_value(material)
    local _struct, custom_fields = expand_point_parameters(params_expr, matval)
    local _checks = quote
        Peridynamics.typecheck_material($(esc(material)))
        Peridynamics.typecheck_params($(esc(material)), $(esc(params)))
    end
    local _pointparam_type = quote
        Peridynamics.point_param_type(::$(esc(material))) = $(esc(params))
    end
    local _material_type = quote
        Peridynamics.material_type(::$(esc(params))) = $(esc(material))
    end
    local _get_pointparams = quote
        function Peridynamics.get_point_params(m::$(esc(material)), p::Dict{Symbol,Any})
            params = extract_parameters(m, $custom_fields, p)
            return $(esc(params))(params...)
        end
    end
    return Expr(:block, _struct, _checks, _pointparam_type, _material_type,
                _get_pointparams)
end

function macrocheck_input_material(material)
    material isa Symbol && return nothing
    (material isa Expr && material.head === :.) && return nothing
    return throw(ArgumentError("argument `$material` is not a valid material input!\n"))
end

function macrocheck_input_params(struct_expr)
    if struct_expr isa Symbol || struct_expr.head != :struct
        msg = "specified input is not a valid point parameter struct expression!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function extract_struct_name(struct_expr::Expr)
    struct_header = struct_expr.args[2]
    struct_header isa Symbol && return struct_header
    struct_header isa Expr || throw(ArgumentError("invalid struct header!\n"))
    if !isa(struct_header, Expr)
        msg = "invalid struct header: $struct_header\n"
        throw(ArgumentError(msg))
    end
    return struct_header.args[1]
end

function extract_material_value(material)
    if material isa Symbol
        return Val(material)
    elseif material isa Expr && material.head === :.
        return Val(material.args[2])
    else
        msg = "invalid material input: $material\n"
        throw(ArgumentError(msg))
    end
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
        msg = "$Param is not a subtype of `AbstractPointParameters`!\n"
        throw(ArgumentError(msg))
    end
    parameters = fieldnames(Param)
    for req_param in required_point_parameter_names(Val(Symbol(Material)))
        if !in(req_param, parameters)
            msg = "required parameter $req_param not found in $(Param)!\n"
            throw(ArgumentError(msg))
        end
    end
    return nothing
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
