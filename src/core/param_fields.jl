#=
This file contains everything that describes *what a point parameter is*, before `@params`
generates a struct from it. It is the parameter counterpart of `core/storage_fields.jl` and
follows the same three ideas:

  * a declaration says as much as possible, so that the macro can generate the struct, the
    constructor, the allowed `material!` keywords and the log methods instead of asking for
    them separately,
  * reuse happens through named, semantic blocks that are included with `@inherit`, so a
    contract is inherited rather than copied,
  * the annotations are never expanded. They are recognized by name while the body is
    parsed, so `@derived` and `Peridynamics.@derived` both work.

It is included next to `core/storage_fields.jl`, whose helpers it reuses, so that parameter
blocks can be declared next to the system, damage model or material family whose contract
they express.
=#

# --------------------------------------------------------------------------------------
# declarations
# --------------------------------------------------------------------------------------

"""
    ParamFieldDecl

$(internal_api_warning())

One parameter declaration of a [`@params`](@ref) or [`@params_fields`](@ref) body.

# Fields

- `name::Symbol`: Name of the parameter, i.e. the field of the generated struct.
- `type::Any`: Declared type of the parameter: a `Type`, [`SimFloat`](@ref) if the
    parameter is declared without a type and follows the float type of the simulation, or a
    type expression over `FT`, e.g. `SArray{NTuple{4,3},FT,4,81}`, which follows it too.
- `kwarg::Symbol`: `material!` keyword the parameter is read from, or `:none` if it is not
    read from a keyword of its own.
- `default::Any`: Expression of the default value, or `nothing` if the keyword is required.
- `provider::Any`: Expression of the call that supplies the parameter, or `nothing` if the
    parameter is not a member of a [`@derived`](@ref) group.
- `source::String`: The right-hand side as it was written, for [`block_table`](@ref) and the
    error messages. It is *not* part of `==`, because two declarations mean the same thing
    when they generate the same code, no matter how they were spelled.
- `label::String`: Label of the simulation-log line of the parameter, or `""` if the
    parameter is not logged.
"""
struct ParamFieldDecl
    name::Symbol
    type::Any
    kwarg::Symbol
    default::Any
    provider::Any
    source::String
    label::String
end

function Base.:(==)(a::ParamFieldDecl, b::ParamFieldDecl)
    return a.name === b.name && a.type == b.type && a.kwarg === b.kwarg &&
           a.default == b.default && a.provider == b.provider && a.label == b.label
end

"""
    ParamFieldsSpec

$(internal_api_warning())

The result of parsing a [`@params`](@ref) or [`@params_fields`](@ref) body: the flattened
[`ParamFieldDecl`](@ref)s in declaration order, and every `material!` keyword they consume,
including the keywords a provider call reads that have no parameter of their own.
"""
struct ParamFieldsSpec
    decls::Vector{ParamFieldDecl}
    kwargs::Vector{Symbol}
end

"""
    param_fields_expr(::Type{T})

$(internal_api_warning())

Return the [`ParamFieldsSpec`](@ref) of a point parameter type or of a parameter block
defined with [`@params_fields`](@ref). This is the table [`@inherit`](@ref) reads at macro
expansion time.
"""
function param_fields_expr(::Type{T}) where {T}
    msg = "`$(T)` does not define any point parameters!\n"
    msg *= "  Only point parameters defined with `@params` and parameter blocks defined "
    msg *= "with `@params_fields` can be inherited from.\n"
    return throw(ArgumentError(msg))
end

"""
    ConstitutiveParameters

$(extension_api_note())

Marker for one field of a [`@params`](@ref) definition or a [`@params_fields`](@ref) block,
e.g. `cm_params::ConstitutiveParameters`. The field holds the parameters the constitutive
model of the material declares with [`@cm_params`](@ref), whatever model that is. The
parameter type is resolved per model when the point parameter type is instantiated, and it
is `Nothing` for a model that declares no parameters. A material whose point parameters
carry this marker therefore supports every constitutive model, parameterized or not,
without knowing any of them.

See also [`DamageParameters`](@ref) and the storage-side analogue `ConstitutiveState`.
"""
struct ConstitutiveParameters end

"""
    DamageParameters

$(extension_api_note())

Marker for one field of a [`@params`](@ref) definition or a [`@params_fields`](@ref) block,
e.g. `dmg_params::DamageParameters`. The field holds the parameters the damage model of the
material declares with [`@dmg_params`](@ref), resolved per model exactly like
[`ConstitutiveParameters`](@ref). The standard fracture parameters `Gc` and `εc` live here:
they belong to [`CriticalStretch`](@ref), not to the material.
"""
struct DamageParameters end

# the type parameters of the generated struct that hold the model parameter types; `CMP`
# comes before `DMP`, mirroring `CMS`/`DMS` of the storage
const CM_PARAM_TYPE_PARAM = :CMP
const DMG_PARAM_TYPE_PARAM = :DMP

@inline is_cm_param_decl(d::ParamFieldDecl) = d.type === ConstitutiveParameters
@inline is_dmg_param_decl(d::ParamFieldDecl) = d.type === DamageParameters
@inline is_model_param_decl(d::ParamFieldDecl) = is_cm_param_decl(d) || is_dmg_param_decl(d)

"""
    constitutive_param_type(model, FT)
    damage_param_type(dmgmodel, FT)

$(internal_api_warning())

Return the point parameter type a model declares with [`@cm_params`](@ref) /
[`@dmg_params`](@ref), instantiated with the float type `FT`, or `Nothing` for a model
without parameters. This is how `point_param_type` resolves the marker fields of a
[`@params`](@ref) definition.
"""
function constitutive_param_type end
constitutive_param_type(model, ::Type) = Nothing
damage_param_type(dmgmodel, ::Type) = Nothing

@doc (@doc constitutive_param_type) damage_param_type

"""
    get_cm_params(model, FT, mat, mat_params, p)
    get_dmg_params(dmgmodel, FT, mat, mat_params, p)

$(internal_api_warning())

Construct the model-owned parameters of one point parameter set, or `nothing` for a model
without parameters. Called by the generated point parameter constructor at the position of
the marker field, with `mat` the material, `mat_params` the `NamedTuple` of every material
parameter declared above the marker and `p` the keyword dictionary of `material!`. The
methods are generated by [`@cm_params`](@ref) / [`@dmg_params`](@ref).
"""
function get_cm_params end
get_cm_params(model, ::Type, mat, mat_params::NamedTuple, p::Dict{Symbol,Any}) = nothing
get_dmg_params(dmgmodel, ::Type, mat, mat_params::NamedTuple, p::Dict{Symbol,Any}) = nothing

@doc (@doc get_cm_params) get_dmg_params

"""
    constitutive_param_kwargs(model)
    damage_param_kwargs(dmgmodel)

$(internal_api_warning())

The `material!` keywords the parameters of a model consume, `()` for a model without
parameters. `material!` accepts the union of the material's own keywords and the keywords
of its models, so a keyword is accepted exactly when something reads it.
"""
function constitutive_param_kwargs end
constitutive_param_kwargs(model) = ()
damage_param_kwargs(dmgmodel) = ()

@doc (@doc constitutive_param_kwargs) damage_param_kwargs

# converting nested model parameters to another float type; a model without parameters has
# `nothing`, a non-generic parameter type converts to itself
convert_nested_params(::Type, ::Nothing) = nothing

function required_model_param(mat_params::NamedTuple, name::Symbol, model)
    haskey(mat_params, name) && return mat_params[name]
    msg = "the parameters of `$(nameof(typeof(model)))` read `$(name)`, which the point "
    msg *= "parameters of the material do not provide!\n"
    msg *= "  A model parameter block sees the material parameters declared above the "
    msg *= "marker field\n"
    msg *= "  (`cm_params::ConstitutiveParameters` / `dmg_params::DamageParameters`). "
    msg *= "Move the marker\n"
    msg *= "  below the declaration of `$(name)`, or use a material whose parameters "
    msg *= "declare it.\n"
    return throw(ArgumentError(msg))
end

@inline is_provided(d::ParamFieldDecl) = !isnothing(d.provider)
@inline is_derived(d::ParamFieldDecl) = isnothing(d.provider) && d.kwarg === :none

@inline function param_field_type(d::ParamFieldDecl)
    is_cm_param_decl(d) && return CM_PARAM_TYPE_PARAM
    is_dmg_param_decl(d) && return DMG_PARAM_TYPE_PARAM
    d.type === SimFloat && return FLOAT_TYPE_PARAM
    return float_type_expr(d.type, FLOAT_TYPE_PARAM)
end

# the name a declaration uses for the float type of the simulation; it is the name of the
# generated type parameter, so `FT` in `SArray{NTuple{4,3},FT,4,81}` is that parameter
const SIM_FLOAT_NAME = :FT

# a parameter follows the float type of the simulation if it is declared without a type or
# with a type expression that names `FT`
@inline function follows_sim_float(d::ParamFieldDecl)
    return d.type === SimFloat || is_float_type_expr(d.type)
end

@inline function is_float_type_expr(type)
    return type isa Expr && in(SIM_FLOAT_NAME, referenced_names(type))
end

@inline function params_use_sim_float(spec::ParamFieldsSpec)
    return any(follows_sim_float, spec.decls)
end

# the type expression with the float type replaced by `ft`: the type parameter of the
# struct and the constructor, or `FT_TO` of the converting constructor
float_type_expr(type, ft) = type
float_type_expr(type::Symbol, ft) = type === SIM_FLOAT_NAME ? ft : type
function float_type_expr(type::Expr, ft)
    return Expr(type.head, (float_type_expr(a, ft) for a in type.args)...)
end

# --------------------------------------------------------------------------------------
# annotations
# --------------------------------------------------------------------------------------

"""
    @kwarg keyword parameter

$(extension_api_note())

Declare the `material!` keyword of a parameter when it differs from the name of the
parameter itself, e.g. `@kwarg gamma_c gammac = 1e-10`.
"""
macro kwarg(keyword, parameter)
    return macro_outside_params_error("@kwarg")
end

"""
    @derived parameter = expression
    @derived (; parameter, ...) = call(args...; keyword, ...)

$(extension_api_note())

Declare a parameter that is computed and is **not** a `material!` keyword of its own. The
second form declares a group of parameters that one call supplies together, such as the six
elastic parameters, which are resolved from any two of six keywords. The call has to return
a `NamedTuple` with one entry per destructured name.

Every right-hand side follows one rule:

| position | meaning |
|:---|:---|
| before `;` | the parameters declared above, `mat`, and anything in scope of the module |
| after `;` | names of `material!` keywords |

```julia
@derived bc = 18 * K / (π * δ^4)
@derived (; δ, rho) = get_discretization_params(; horizon, rho)
```

The keywords written after `;` are exactly the `material!` keywords the block consumes, so
`allowed_material_kwargs` cannot disagree with the call that reads them. A keyword the user
did not pass is not forwarded, which is why the provider decides on its own whether a
keyword is required:

```julia
get_discretization_params(; horizon, rho)        # both required
get_elastic_params(; E=nothing, nu=nothing, …)   # any two of six
```

A body may declare any number of groups, in any order. Without this annotation a parameter
with a default expression is a keyword that defaults to that expression.
"""
macro derived(parameter)
    return macro_outside_params_error("@derived")
end

"""
    @log "label" parameter
    @log "label" declaration

$(extension_api_note())

Write a parameter to the simulation log under `label`. The first form labels a parameter
that is already declared above, which is how a member of a [`@derived`](@ref) group is
logged, and the second declares and labels in one line. Without this annotation the parameter
does not appear in the log, which is the default of `log_param_property`.

```julia
@log "yield stress" sigma_y = Inf
@log "bond constant" @derived bc = 18 * K / (π * δ^4)

@derived (; δ, rho) = get_discretization_params(; horizon, rho)
@log "horizon" δ
```
"""
macro log(label, parameter)
    return macro_outside_params_error("@log")
end

"""
    material_kwargs(p, names)

$(internal_api_warning())

Return the `material!` keywords `names` that are present in `p` as a `NamedTuple`. This is
what a provider call of a [`@derived`](@ref) declaration is passed: a keyword the user did
not give is not forwarded, so the provider decides whether it is required.
"""
function material_kwargs(p::Dict{Symbol,Any}, names::NTuple{N,Symbol}) where {N}
    present = filter(name -> haskey(p, name), names)
    return NamedTuple{present}(map(name -> p[name], present))
end

function macro_outside_params_error(name::AbstractString)
    msg = "`$(name)` is only allowed inside a `@params` or `@params_fields` definition!\n"
    return :(throw(ArgumentError($msg)))
end

# --------------------------------------------------------------------------------------
# `@params_fields`
# --------------------------------------------------------------------------------------

"""
    @params_fields name begin ... end

$(extension_api_note())

Define a reusable block of point parameter declarations that can be included into a
[`@params`](@ref) definition with [`@inherit`](@ref). The body accepts exactly the same
declarations as `@params`, including `@inherit` itself.

The macro defines a marker type `name` and registers the flattened declarations with
[`param_fields_expr`](@ref). A block emits no method of its own, in particular no
`log_param_property`, because several blocks may describe the same parameter, e.g.
`ElasticParameters` and `BBElasticParameters`. The `@log` labels travel with the
declarations and are emitted by the [`@params`](@ref) definition that defines the type.

# Example
```julia
Peridynamics.@params_fields DiscretizationParameters begin
    @derived (; δ, rho) = get_discretization_params(; horizon, rho)
    @log "horizon" δ
    @log "density" rho
end
```

The blocks this package ships are [`DiscretizationParameters`](@ref),
[`ElasticParameters`](@ref), [`BBElasticParameters`](@ref),
[`BondHorizonParameters`](@ref), [`InteractionParameters`](@ref) and
[`StandardParameters`](@ref), each documented with the parameters it exposes;
[`FractureParameters`](@ref) belongs to the damage model and is inherited inside a
[`@dmg_params`](@ref) declaration. The point parameters of a material defined with
[`@params`](@ref) can be inherited the same way, e.g. `@inherit BBPointParameters`.
"""
macro params_fields(name, block)
    macrocheck_input_params_fields_name(name)
    macrocheck_input_params_fields_block(block)
    local _spec = get_param_decls(block.args, __module__)
    local _struct = quote
        struct $(esc(name)) <: Peridynamics.AbstractPointParameterFields end
    end
    local _fields_expr = quote
        function Peridynamics.param_fields_expr(::Base.Type{$(esc(name))})
            return $(QuoteNode(_spec))
        end
    end
    # the docstring is attached last, so that a `$(block_table(Name))` in it can already read
    # the declarations that were just registered
    local _doc = quote
        Base.@__doc__ $(esc(name))
    end
    return Expr(:block, _struct, _fields_expr, _doc)
end

function macrocheck_input_params_fields_name(name)
    name isa Symbol && return nothing
    msg = "argument `$name` is not a valid parameter block name!\n"
    return throw(ArgumentError(msg))
end

function macrocheck_input_params_fields_block(block)
    (block isa Expr && block.head === :block) && return nothing
    msg = "specified input is not a valid parameter block expression!\n"
    return throw(ArgumentError(msg))
end

# --------------------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------------------

"""
    get_param_decls(block_args, mod)

$(internal_api_warning())

Parse the body of a [`@params`](@ref) or [`@params_fields`](@ref) definition into a
[`ParamFieldsSpec`](@ref), expanding every [`@inherit`](@ref) and applying the merge rules
documented there.

Types are resolved in `mod` while the macro is expanded, so an inherited declaration always
keeps the meaning it has in the module that declared it. Default and provider expressions are
escaped instead, because they are evaluated in the generated constructor, where they see the
material `mat`, the keyword dictionary `p` and the parameters declared before them.
"""
function get_param_decls(block_args, mod::Module)
    decls = Vector{ParamFieldDecl}()
    kwargs = Vector{Symbol}()
    sources = Dict{Symbol,String}()
    for field in block_args
        field isa LineNumberNode && continue
        if !isa(field, Expr) && !isa(field, Symbol)
            throw(ArgumentError("unexpected parameter declaration: $field\n"))
        end
        macro_name = (field isa Expr && field.head === :macrocall) ?
                     get_macro_name(field.args[1]) : nothing
        if macro_name === Symbol("@inherit")
            for inherit_expr in get_inherit_arguments(field)
                source = string(inherit_expr)
                inherited = resolve_inherited_params(mod, inherit_expr)
                for decl in inherited.decls
                    add_param_decl!(decls, sources, decl, source, true)
                end
                for kwarg in inherited.kwargs
                    in(kwarg, kwargs) || push!(kwargs, kwarg)
                end
            end
        elseif attach_label!(decls, field)
            continue
        else
            for decl in parse_param_decls(field, mod, kwargs)
                add_param_decl!(decls, sources, decl, "the definition itself", false)
                decl.kwarg === :none || in(decl.kwarg, kwargs) || push!(kwargs, decl.kwarg)
            end
        end
    end
    check_param_order(decls, mod)
    check_single_model_params(decls)
    return ParamFieldsSpec(decls, kwargs)
end

#=
`@log "label" name` with a bare name labels a parameter that is already declared, which is
how a member of a `@derived` group is logged. The same syntax on a name that is not declared
yet declares a required keyword and labels it, so the two forms never overlap.
=#
function attach_label!(decls::Vector{ParamFieldDecl}, field)
    field isa Expr && field.head === :macrocall || return false
    get_macro_name(field.args[1]) === Symbol("@log") || return false
    args = filter(x -> !isa(x, LineNumberNode), field.args[2:end])
    (length(args) == 2 && args[1] isa AbstractString && args[2] isa Symbol) || return false
    idx = findfirst(d -> d.name === args[2], decls)
    isnothing(idx) && return false
    d = decls[idx]
    decls[idx] = ParamFieldDecl(d.name, d.type, d.kwarg, d.default, d.provider, d.source,
                                String(args[1]))
    return true
end

function resolve_inherited_params(mod::Module, name_expr)
    inherited_type = resolve_in_mod_or_peridynamics(mod, name_expr)
    if isnothing(inherited_type)
        msg = "cannot resolve `$(name_expr)` of `@inherit`!\n"
        msg *= "  It has to be defined by an earlier top-level statement and it has to be "
        msg *= "resolvable in `$(mod)` or in `Peridynamics`.\n"
        throw(ArgumentError(msg))
    end
    if !isa(inherited_type, Type)
        throw(ArgumentError("`$(name_expr)` of `@inherit` is not a type!\n"))
    end
    # `invokelatest` because the method of `param_fields_expr` for `inherited_type` was
    # defined by a top-level statement that is newer than this macro
    return Base.invokelatest(param_fields_expr, inherited_type)
end

function parse_param_decls(field, mod::Module, kwargs::Vector{Symbol})
    (; expr, kwarg, derived, label) = peel_param_annotations(field)
    is_param_group(expr) || return [parse_param_decl(expr, kwarg, derived, label, mod,
                                                     kwargs)]
    derived || throw(ArgumentError(group_needs_derived_msg(expr)))
    isnothing(kwarg) || throw(ArgumentError(group_annotation_msg(expr, "`@kwarg`")))
    isempty(label) || throw(ArgumentError(group_annotation_msg(expr, "a `@log` label")))
    return parse_param_group(expr, mod, kwargs)
end

function parse_param_decl(expr, kwarg, derived, label, mod::Module, kwargs::Vector{Symbol})
    (; name, type, default) = parse_param_field(expr, mod)
    if type === ConstitutiveParameters || type === DamageParameters
        check_marker_decl(expr, type, kwarg, derived, label, default)
        return ParamFieldDecl(name, type, :none, nothing, nothing, "", "")
    end
    param_kwarg = derived ? :none : (isnothing(kwarg) ? name : kwarg)
    isnothing(default) && return ParamFieldDecl(name, type, param_kwarg, nothing, nothing,
                                                "", label)
    check_no_keyword_dict(default)
    esc_default = resolve_param_expr(mod, read_material_kwargs(default, kwargs))
    return ParamFieldDecl(name, type, param_kwarg, esc_default, nothing, string(default),
                          label)
end

#=
A group is `(; a, b) = call(...)`: one call supplies several parameters at once. Every member
shares the same resolved provider expression, so `param_constructor_body` emits one call for
the whole group.
=#
@inline function is_param_group(expr)
    expr isa Expr && expr.head === :(=) || return false
    lhs = expr.args[1]
    return lhs isa Expr && lhs.head === :tuple && length(lhs.args) == 1 &&
           lhs.args[1] isa Expr && lhs.args[1].head === :parameters
end

function parse_param_group(expr, mod::Module, kwargs::Vector{Symbol})
    call = expr.args[2]
    if !(call isa Expr && call.head === :call)
        throw(ArgumentError(group_needs_call_msg(expr)))
    end
    check_no_keyword_dict(call)
    provider = resolve_param_expr(mod, read_material_kwargs(call, kwargs))
    source = string(call)
    decls = Vector{ParamFieldDecl}()
    for member in expr.args[1].args[1].args
        if member isa Symbol
            check_param_name(member)
            push!(decls, ParamFieldDecl(member, SimFloat, :none, nothing, provider, source,
                                        ""))
        elseif member isa Expr && member.head === :(::) && member.args[1] isa Symbol
            check_param_name(member.args[1])
            type = resolve_param_type(mod, member.args[2])
            push!(decls, ParamFieldDecl(member.args[1], type, :none, nothing, provider,
                                        source, ""))
        else
            throw(ArgumentError(group_member_msg(member)))
        end
    end
    isempty(decls) && throw(ArgumentError(group_needs_call_msg(expr)))
    return decls
end

function group_needs_derived_msg(expr)
    msg = "a group of parameters has to be declared with `@derived`, got: $(expr)\n"
    msg *= "  Every parameter a call supplies is computed and is not a `material!` keyword "
    msg *= "of its own, e.g.\n"
    msg *= "        @derived (; δ, rho) = get_discretization_params(; horizon, rho)\n"
    return msg
end

function group_needs_call_msg(expr)
    msg = "the right-hand side of a parameter group is not a call, got: $(expr)\n"
    msg *= "  A group is supplied by one call that returns a `NamedTuple` with one entry "
    msg *= "per destructured name.\n"
    return msg
end

function group_member_msg(member)
    msg = "unexpected member of a parameter group: $(member)\n"
    msg *= "  Members are declared as `name` or `name::Type`.\n"
    return msg
end

function group_annotation_msg(expr, what)
    msg = "a parameter group cannot declare $(what), got: $(expr)\n"
    msg *= "  It applies to a single parameter. Annotate the members individually, e.g. "
    msg *= "`@log \"critical stretch\" εc` below the group.\n"
    return msg
end

#=
The `material!` keywords a call reads are its keyword arguments, written in shorthand. They
are collected here and replaced by the keywords that are actually present, so a provider with
a required keyword argument reports a missing one itself, with an `UndefKeywordError`. An
explicit `key = value` is rejected: it would make "after `;` is a `material!` keyword" a rule
with an exception, and a parameter of an earlier declaration belongs in a positional argument.
=#
function read_material_kwargs(expr, kwargs::Vector{Symbol})
    expr isa Expr || return expr
    expr.head === :call || return Expr(expr.head,
                                       (read_material_kwargs(a, kwargs)
                                        for a in expr.args)...)
    args = Any[expr.args[1]]
    for (i, arg) in enumerate(expr.args[2:end])
        if arg isa Expr && arg.head === :parameters && i == 1
            push!(args, material_kwargs_expr(arg, expr, kwargs))
        elseif arg isa Expr && arg.head === :kw
            throw(ArgumentError(explicit_kwarg_msg(expr, arg)))
        else
            push!(args, read_material_kwargs(arg, kwargs))
        end
    end
    return Expr(:call, args...)
end

function material_kwargs_expr(params::Expr, call::Expr, kwargs::Vector{Symbol})
    names = Symbol[]
    for arg in params.args
        arg isa Symbol || throw(ArgumentError(explicit_kwarg_msg(call, arg)))
        push!(names, arg)
        in(arg, kwargs) || push!(kwargs, arg)
    end
    selection = Expr(:call, material_kwargs, :p, Expr(:tuple, QuoteNode.(names)...))
    return Expr(:parameters, Expr(:..., selection))
end

function explicit_kwarg_msg(call, arg)
    shown = arg isa Expr && arg.head === :kw ? "$(arg.args[1]) = $(arg.args[2])" : "$(arg)"
    msg = "unexpected keyword argument `$(shown)` in: $(call)\n"
    msg *= "  The keyword arguments of a call are the names of the `material!` keywords it "
    msg *= "reads, so they are written in shorthand:\n"
    msg *= "        get_bond_horizon(δ; bond_horizon)\n"
    msg *= "  Pass the parameters declared above positionally.\n"
    return msg
end

#=
`p`, the keyword dictionary of `material!`, used to be spliced into every declaration. The
keywords are named at the call that reads them now, so the name has no meaning any more and
is rejected rather than silently resolving to something else.
=#
function check_no_keyword_dict(expr)
    in(:p, referenced_names(expr)) || return nothing
    msg = "`p` is not available in a parameter declaration: $(expr)\n"
    msg *= "  Name the `material!` keywords a call reads as its keyword arguments:\n"
    msg *= "        @derived (; δb) = get_bond_horizon(δ; bond_horizon)\n"
    return throw(ArgumentError(msg))
end

function peel_param_annotations(field)
    kwarg = nothing
    derived = false
    label = ""
    expr = field
    while expr isa Expr && expr.head === :macrocall
        macro_name = get_macro_name(expr.args[1])
        args = filter(x -> !isa(x, LineNumberNode), expr.args[2:end])
        if macro_name === Symbol("@derived")
            length(args) == 1 || throw(ArgumentError(annotation_msg("@derived", expr)))
            derived = true
            expr = only(args)
        elseif macro_name === Symbol("@kwarg")
            (length(args) == 2 && args[1] isa Symbol) ||
                throw(ArgumentError(annotation_msg("@kwarg", expr)))
            kwarg = args[1]
            expr = args[2]
        elseif macro_name === Symbol("@log")
            (length(args) == 2 && args[1] isa AbstractString) ||
                throw(ArgumentError(annotation_msg("@log", expr)))
            label = String(args[1])
            expr = args[2]
        else
            msg = "unknown annotation in a parameter declaration: $(expr)\n"
            msg *= "  Known annotations are @derived, @kwarg, @log and @inherit.\n"
            throw(ArgumentError(msg))
        end
    end
    return (; expr, kwarg, derived, label)
end

function annotation_msg(name, expr)
    msg = "malformed `$(name)` annotation: $(expr)\n"
    name == "@derived" && (msg *= "  Expected `@derived name = expression`.\n")
    name == "@kwarg" && (msg *= "  Expected `@kwarg keyword name`.\n")
    name == "@log" && (msg *= "  Expected `@log \"label\" name`.\n")
    return msg
end

#=
A marker is one bare field, `cm_params::ConstitutiveParameters` or
`dmg_params::DamageParameters`: the model owns the declarations behind it, so nothing of
the declaration language applies to the field itself.
=#
function check_marker_decl(expr, type, kwarg, derived, label, default)
    problem = if !isnothing(default)
        "cannot have a default value"
    elseif derived
        "cannot be `@derived`"
    elseif !isnothing(kwarg)
        "cannot have a `@kwarg` keyword"
    elseif !isempty(label)
        "cannot carry a `@log` label"
    else
        nothing
    end
    isnothing(problem) && return nothing
    msg = "the marker field `$(expr)` $(problem)!\n"
    msg *= "  The $(marker_model_name(type)) owns the declarations behind the marker, so "
    msg *= "the marker itself\n  declares nothing.\n"
    return throw(ArgumentError(msg))
end

marker_model_name(::Type{ConstitutiveParameters}) = "constitutive model"
marker_model_name(::Type{DamageParameters}) = "damage model"

function check_single_model_params(decls::Vector{ParamFieldDecl})
    for predicate in (is_cm_param_decl, is_dmg_param_decl)
        found = filter(predicate, decls)
        length(found) <= 1 && continue
        names = join(("`" * string(d.name) * "`" for d in found), ", ")
        model = marker_model_name(found[1].type)
        msg = "point parameters can carry the parameters of the $(model) only once, "
        msg *= "found $(names)!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function parse_param_field(expr, mod::Module)
    default = nothing
    decl = expr
    if decl isa Expr && decl.head === :(=)
        default = decl.args[2]
        decl = decl.args[1]
    end
    if decl isa Symbol
        check_param_name(decl)
        # an untyped parameter follows the float type of the simulation
        return (; name=decl, type=SimFloat, default)
    end
    if decl isa Expr && decl.head === :(::) && decl.args[1] isa Symbol
        check_param_name(decl.args[1])
        return (; name=decl.args[1], type=resolve_param_type(mod, decl.args[2]), default)
    end
    return throw(ArgumentError("unexpected parameter declaration: $expr\n"))
end

# `FT` names the float type of the simulation inside a definition, so it cannot name a
# parameter
function check_param_name(name::Symbol)
    name === SIM_FLOAT_NAME || return nothing
    msg = "`$(name)` cannot be the name of a parameter!\n"
    msg *= "  Inside a parameter definition `$(SIM_FLOAT_NAME)` stands for the float type "
    msg *= "of the simulation.\n"
    return throw(ArgumentError(msg))
end

#=
A provider call and a default expression are evaluated in the generated constructor, so they
are escaped and see the parameters, `mat` and `p` of that constructor. The *functions* they
call have to keep the meaning they have in the module of the declaration, though, because an
inherited declaration is spliced into a definition in another module, where the name of a
function of this package does not resolve. So every callee is resolved while the macro is
expanded and the resolved function is interpolated, the same way `@storage` interpolates the
resolved type of a field. Everything that is not a callee stays untouched, because it may be
a parameter of the generated constructor.
=#
function resolve_param_expr(mod::Module, expr)
    return esc(resolve_callees(mod, expr))
end

resolve_callees(::Module, x) = x

function resolve_callees(mod::Module, expr::Expr)
    args = copy(expr.args)
    start = 1
    if expr.head === :call && !isempty(args)
        args[1] = resolve_callee(mod, args[1])
        start = 2
    end
    for i in start:length(args)
        args[i] = resolve_callees(mod, args[i])
    end
    return Expr(expr.head, args...)
end

function resolve_callee(mod::Module, callee)
    resolved = resolve_in_mod_or_peridynamics(mod, callee)
    isa(resolved, Function) && return resolved
    # A provider may well be defined by a *later* top-level statement of the same module,
    # e.g. a block declared next to the concept it describes and the function below it. Such
    # a name cannot be resolved while the macro is expanded, so it is qualified with the
    # module of the declaration instead and resolves when the generated constructor runs.
    callee isa Symbol && return Expr(:., mod, QuoteNode(callee))
    return callee
end

function resolve_param_type(mod::Module, type_expr)
    is_float_type_expr(type_expr) && return resolve_float_type_expr(mod, type_expr)
    resolved = resolve_in_mod_or_peridynamics(mod, type_expr)
    isa(resolved, Type) && return resolved
    msg = "cannot resolve the type `$(type_expr)` of a point parameter!\n"
    msg *= "  It has to be resolvable in `$(mod)` or in `Peridynamics`. Declare the "
    msg *= "parameter without a type to let it follow the float type of the simulation, "
    msg *= "or\n  build the type from `$(SIM_FLOAT_NAME)`, e.g. "
    msg *= "`SArray{NTuple{4,3},$(SIM_FLOAT_NAME),4,81}`.\n"
    return throw(ArgumentError(msg))
end

#=
A type expression that names `FT` cannot be evaluated while the macro is expanded, because
`FT` is the type parameter of the struct that is about to be generated. So it stays an
expression, and every other name in it is resolved to the module that defines it, the way a
provider call is resolved: an inherited declaration keeps its meaning in another module, and
the expression evaluates wherever the generated struct and constructors are.
=#
function resolve_float_type_expr(mod::Module, type_expr)
    return resolve_type_names(mod, type_expr, type_expr)
end

resolve_type_names(::Module, x, type_expr) = x

function resolve_type_names(mod::Module, name::Symbol, type_expr)
    name === SIM_FLOAT_NAME && return name
    isnothing(try_eval(mod, name)) || return GlobalRef(mod, name)
    isnothing(try_eval(Peridynamics, name)) || return GlobalRef(Peridynamics, name)
    msg = "cannot resolve `$(name)` in the type `$(type_expr)` of a point parameter!\n"
    msg *= "  Every name but `$(SIM_FLOAT_NAME)` has to be resolvable in `$(mod)` or in "
    msg *= "`Peridynamics`.\n"
    return throw(ArgumentError(msg))
end

function resolve_type_names(mod::Module, expr::Expr, type_expr)
    # `a.b` resolves `a`, never `b`, and a quoted name is data
    expr.head === :quote && return expr
    if expr.head === :.
        return Expr(:., resolve_type_names(mod, expr.args[1], type_expr),
                    expr.args[2:end]...)
    end
    return Expr(expr.head, (resolve_type_names(mod, a, type_expr) for a in expr.args)...)
end

# the type as it was written, for the tables and the error messages
type_expr_string(type) = type_msg(type)
type_expr_string(type::Expr) = string(unresolve_type_names(type))

unresolve_type_names(x) = x
unresolve_type_names(ref::GlobalRef) = ref.name
function unresolve_type_names(expr::Expr)
    return Expr(expr.head, (unresolve_type_names(a) for a in expr.args)...)
end

#=
A right-hand side sees the parameters declared above it, so a reference to one declared below
is the mistake this rule makes easy: it is the reason a block has to follow the block that
provides what it reads. The names of the calls have been resolved to function objects and the
`material!` keywords to a selection call by now, so a bare symbol left in the expression is a
value. Names that resolve in the module are globals and are left alone.
=#
function check_param_order(decls::Vector{ParamFieldDecl}, mod::Module)
    for (i, decl) in enumerate(decls)
        exprs = (decl.default, decl.provider)
        for expr in exprs
            isnothing(expr) && continue
            for name in referenced_names(expr)
                idx = findfirst(d -> d.name === name, decls)
                (isnothing(idx) || idx <= i) && continue
                isnothing(resolve_in_mod_or_peridynamics(mod, name)) || continue
                throw(ArgumentError(param_order_msg(decl, decls[idx])))
            end
        end
    end
    return nothing
end

function param_order_msg(decl::ParamFieldDecl, referenced::ParamFieldDecl)
    msg = "the parameter `$(decl.name)` reads `$(referenced.name)`, which is declared "
    msg *= "below it!\n"
    msg *= "  A right-hand side sees only the parameters declared above it. Move the "
    msg *= "declaration of `$(referenced.name)`, or the `@inherit` that contributes it, "
    msg *= "above `$(decl.name)`.\n"
    return msg
end

referenced_names(x) = collect_referenced_names!(Set{Symbol}(), x)

collect_referenced_names!(names::Set{Symbol}, x) = names
collect_referenced_names!(names::Set{Symbol}, x::Symbol) = push!(names, x)

function collect_referenced_names!(names::Set{Symbol}, expr::Expr)
    # `a.b` references `a`, never `b`, and a quoted name is data
    expr.head === :quote && return names
    args = expr.head === :. ? expr.args[1:1] : expr.args
    for arg in args
        collect_referenced_names!(names, arg)
    end
    return names
end

function add_param_decl!(decls::Vector{ParamFieldDecl}, sources::Dict{Symbol,String},
                         decl::ParamFieldDecl, source::AbstractString, from_inherit::Bool)
    idx = findfirst(d -> d.name === decl.name, decls)
    if isnothing(idx)
        push!(decls, decl)
        sources[decl.name] = source
        return nothing
    end
    if !from_inherit
        # a parameter of the definition itself overrides an inherited one in place
        decls[idx] = decl
        sources[decl.name] = source
        return nothing
    end
    decls[idx] == decl && return nothing
    msg = "conflicting declarations of the point parameter `$(decl.name)`!\n"
    msg *= "  inherited from $(sources[decl.name]): $(param_decl_msg(decls[idx]))\n"
    msg *= "  inherited from $(source): $(param_decl_msg(decl))\n"
    msg *= "  Declare the parameter in the definition itself to resolve the conflict.\n"
    return throw(ArgumentError(msg))
end

function param_decl_msg(d::ParamFieldDecl)
    msg = isempty(d.label) ? "" : "@log \"$(d.label)\" "
    is_derived(d) && (msg *= "@derived ")
    (d.kwarg !== :none && d.kwarg !== d.name) && (msg *= "@kwarg $(d.kwarg) ")
    msg *= string(d.name)
    d.type === SimFloat || (msg *= "::$(type_expr_string(d.type))")
    isnothing(d.default) || (msg *= " = $(d.source)")
    is_provided(d) && (msg *= " (from $(d.source))")
    return msg
end

# --------------------------------------------------------------------------------------
# code generation shared by `@params_fields` and `@params`
# --------------------------------------------------------------------------------------

#=
The generated constructor is the constructor that used to be written by hand: one
call per `@derived` group, one keyword lookup per keyword parameter, one expression per
derived parameter, then the positional call. The parameter variables are escaped, so that the
escaped default and provider expressions of the definition see them, which is what lets a
default read the parameters declared before it, and what lets a provider call read `mat` and
the selected `material!` keywords.
=#
function param_constructor_body(spec::ParamFieldsSpec, sim_float::Bool)
    body = Any[]
    provider, provider_var = nothing, nothing
    names_above = Symbol[]
    ft_expr = sim_float ? FLOAT_TYPE_PARAM : :(Peridynamics.default_float_type())
    for decl in spec.decls
        if is_model_param_decl(decl)
            # the model constructs its own parameters, so this is the parameter-side
            # `init_constitutive_state`: it sees every material parameter declared above
            # the marker as a `NamedTuple`, and the keyword dictionary of `material!`
            getter = is_cm_param_decl(decl) ? :(Peridynamics.get_cm_params) :
                     :(Peridynamics.get_dmg_params)
            model = is_cm_param_decl(decl) ?
                    :(Peridynamics.get_constitutive_model($(esc(:mat)))) :
                    :(Peridynamics.get_dmgmodel($(esc(:mat))))
            nt = Expr(:tuple, Expr(:parameters, (esc(n) for n in names_above)...))
            value = Expr(:call, getter, model, ft_expr, esc(:mat), nt, esc(:p))
            push!(body, Expr(:(=), esc(decl.name), value))
            continue
        end
        if is_provided(decl)
            if isnothing(provider) || provider != decl.provider
                provider = decl.provider
                provider_var = Symbol("from_", decl.name)
                push!(body, Expr(:(=), provider_var, provider))
            end
            value = :(Base.getproperty($(provider_var), $(QuoteNode(decl.name))))
        elseif is_derived(decl)
            value = decl.default
        elseif isnothing(decl.default)
            value = :(Peridynamics.required_param($(esc(:p)), $(QuoteNode(decl.kwarg))))
        else
            value = Expr(:if, :(Base.haskey($(esc(:p)), $(QuoteNode(decl.kwarg)))),
                         :($(esc(:p))[$(QuoteNode(decl.kwarg))]), decl.default)
        end
        converted = :(Base.convert($(param_field_type(decl)), $(value)))
        push!(body, Expr(:(=), esc(decl.name), converted))
        push!(names_above, decl.name)
    end
    return body
end

"""
    required_param(p, kwarg)

$(internal_api_warning())

Return the value of a required `material!` keyword, or throw an `UndefKeywordError` naming
it. This is how a point parameter that is declared without a default value is read.
"""
function required_param(p::Dict{Symbol,Any}, kwarg::Symbol)
    haskey(p, kwarg) || throw(UndefKeywordError(kwarg))
    return p[kwarg]
end

function log_material_parameters(param::P; indentation::Int=2) where {P}
    msg = ""
    for key in fieldnames(P)
        value = Base.getfield(param, key)
        if value isa Union{AbstractConstitutiveParameters,AbstractDamageParameters}
            # a marker field: the model's parameters log themselves, with the labels the
            # model declared
            msg *= log_material_parameters(value; indentation)
        else
            msg *= log_param_property(Val(key), param; indentation)
        end
    end
    return msg
end

# a parameter without a `@log` label does not appear in the simulation log, which is why
# this is the fallback and not an error
function log_param_property(::Val{S}, param; indentation) where {S}
    return ""
end

#=
A label belongs to the point parameter *type*, not to a block and not to a material, so the
generated method dispatches on the type and is emitted only by the definition that defines
it. A block emits nothing, because several blocks describe the same parameter with the same
label, e.g. `ElasticParameters` and `BBElasticParameters`, and emitting from both would
redefine the same method. A definition that only adds a constructor emits nothing either,
because the type already answers.
=#
function param_log_methods(spec::ParamFieldsSpec, params_type)
    return [
        quote
            function Peridynamics.log_param_property(::Base.Val{$(QuoteNode(decl.name))},
                                                     param::$(esc(params_type));
                                                     indentation)
                return Peridynamics.msg_qty($(decl.label),
                                            Base.getfield(param, $(QuoteNode(decl.name)));
                                            indentation)
            end
        end
        for decl in spec.decls if !isempty(decl.label)
    ]
end

# --------------------------------------------------------------------------------------
# flat reads over the marker fields
# --------------------------------------------------------------------------------------

#=
`params.Gc` reads like a field even though `Gc` lives in the parameters of the damage
model: `getproperty` forwards a name that is not a field into the marker fields. The
`@generated` function computes the exact `getfield` chain per (type, name), so a property
literal compiles to a direct load and the forwarding costs nothing, which is asserted by the
"parameter property forwarding" item of `test/perf/perf.jl`. Everything generated by the
macros reads fields with `getfield`, so the forwarding cannot recurse.
=#
@inline function Base.getproperty(params::AbstractPointParameters, name::Symbol)
    return get_param_property(params, Val(name))
end

@generated function get_param_property(params::P,
                                       ::Val{S}) where {P<:AbstractPointParameters,S}
    S in fieldnames(P) && return :(Base.getfield(params, $(QuoteNode(S))))
    exprs, owners = Any[], Symbol[]
    for field in fieldnames(P)
        FT = fieldtype(P, field)
        FT <: Union{AbstractConstitutiveParameters,AbstractDamageParameters} || continue
        S in fieldnames(FT) || continue
        push!(owners, field)
        push!(exprs, :(Base.getfield(Base.getfield(params, $(QuoteNode(field))),
                                     $(QuoteNode(S)))))
    end
    length(owners) == 1 && return exprs[1]
    # not found anywhere: `getfield` throws the native no-field error of this Julia version
    isempty(owners) && return :(Base.getfield(params, $(QuoteNode(S))))
    return :(Peridynamics.ambiguous_param_property($(QuoteNode(S)), params,
                                                   $(Tuple(owners))))
end

function ambiguous_param_property(name::Symbol, @nospecialize(params), owners)
    msg = "the parameter `$(name)` exists in $(join(("`$(o)`" for o in owners), " and "))"
    msg *= " of `$(nameof(typeof(params)))`!\n"
    msg *= "  Read it from the model parameters directly, e.g. "
    msg *= "`params.$(first(owners)).$(name)`.\n"
    return throw(ArgumentError(msg))
end

function Base.propertynames(params::AbstractPointParameters, private::Bool=false)
    return param_property_names(params)
end

# the fields plus everything the marker fields hold, which is what `getproperty` accepts
@generated function param_property_names(params::P) where {P<:AbstractPointParameters}
    names = Symbol[]
    for field in fieldnames(P)
        push!(names, field)
        FT = fieldtype(P, field)
        FT <: Union{AbstractConstitutiveParameters,AbstractDamageParameters} || continue
        append!(names, fieldnames(FT))
    end
    return Tuple(names)
end

# the flat parameter names for display: the marker fields replaced by what they hold
@generated function flat_param_property_names(params::P) where {P<:AbstractPointParameters}
    names = Symbol[]
    for field in fieldnames(P)
        FT = fieldtype(P, field)
        if FT <: Union{AbstractConstitutiveParameters,AbstractDamageParameters}
            append!(names, fieldnames(FT))
        elseif FT === Nothing
            continue
        else
            push!(names, field)
        end
    end
    return Tuple(names)
end

# --------------------------------------------------------------------------------------
# struct and conversion code generation, shared by `@params` and the model macros
# --------------------------------------------------------------------------------------

#=
The struct of a parameter definition: parametric in the float type when any declaration
follows the simulation float, plus one unbounded type parameter per marker field, because a
nested model parameter type is whatever the model declares, and `Nothing` for a model without
parameters. `CMP` comes before `DMP`.
=#
function params_struct_expr(params_expr, name, supertype, spec, sim_float)
    fields = [Expr(:(::), d.name, param_field_type(d)) for d in spec.decls]
    type_params = Any[]
    sim_float && push!(type_params, Expr(:(<:), FLOAT_TYPE_PARAM, :(Base.Real)))
    any(is_cm_param_decl, spec.decls) && push!(type_params, CM_PARAM_TYPE_PARAM)
    any(is_dmg_param_decl, spec.decls) && push!(type_params, DMG_PARAM_TYPE_PARAM)
    header = if isempty(type_params)
        Expr(:(<:), name, supertype)
    else
        Expr(:(<:), Expr(:curly, name, type_params...), supertype)
    end
    struct_expr = Expr(:struct, params_expr.args[1], header, Expr(:block, fields...))
    return quote
        $(struct_expr)
    end
end

# converting a whole set of point parameters to another float type is what a simulation in
# `Float32` needs, and it is one method because every parameter converts on its own; a
# marker field converts through the model parameter type it holds
function params_convert_expr(name, spec, sim_float)
    sim_float || return Expr(:block)
    # with the fully implicit positional call the float type is inferred from the values,
    # so every simulation-float field converts explicitly
    getfields = [if is_model_param_decl(d)
                     :(Peridynamics.convert_nested_params(FT_TO,
                                                          Base.getfield(param,
                                                                        $(QuoteNode(d.name)))))
                 elseif follows_sim_float(d)
                     :(Base.convert($(float_type_expr(param_field_type(d), :FT_TO)),
                                    Base.getfield(param, $(QuoteNode(d.name)))))
                 else
                     :(Base.getfield(param, $(QuoteNode(d.name))))
                 end
                 for d in spec.decls]
    instantiation = Expr(:curly, esc(name), :FT_TO)
    # with markers the positional call is the fully implicit one, inferring the model types
    positional = any(is_model_param_decl, spec.decls) ? esc(name) : instantiation
    return quote
        function $(instantiation)(param::$(esc(name))) where {FT_TO}
            return $(positional)($(getfields...))
        end
    end
end

# whether a `@params` type carries a marker field, so that a setup-time check can name the
# missing marker when a model brings parameters the material has no place for
has_cm_param_marker(::Type) = false
has_dmg_param_marker(::Type) = false

function params_marker_exprs(name, spec)
    exprs = Any[]
    if any(is_cm_param_decl, spec.decls)
        push!(exprs,
              :(Peridynamics.has_cm_param_marker(::Base.Type{<:$(esc(name))}) = true))
    end
    if any(is_dmg_param_decl, spec.decls)
        push!(exprs,
              :(Peridynamics.has_dmg_param_marker(::Base.Type{<:$(esc(name))}) = true))
    end
    return exprs
end

# --------------------------------------------------------------------------------------
# `@cm_params` and `@dmg_params`: the parameters a model owns
# --------------------------------------------------------------------------------------

"""
    @cm_params model struct ModelParameters ... end

$(extension_api_note())

Declare the point parameters a constitutive model owns, with the declaration language of
[`@params_fields`](@ref). The macro generates the parameter struct (parametric in the float
type, `isbits`), its constructor, the registration that resolves a
`cm_params::ConstitutiveParameters` marker field of a material to this type, the
`material!` keywords the declarations consume, and the simulation-log methods of the
`@log` labels. So a model registers its own keywords without the material knowing them.

```julia
struct LinearHardening <: Peridynamics.AbstractConstitutiveModel end

Peridynamics.@cm_params LinearHardening struct LinearHardeningParameters
    @log "yield stress" sigma_y
    @log "hardening modulus" H = 0.0
end
```

Inside the block the model instance is available as `model` and the material as `mat`,
and a declaration can read every material parameter declared above the marker field, e.g.
the horizon `δ` or the bulk modulus `K`. A material whose point parameters carry the marker
`cm_params::ConstitutiveParameters` supports the model; its parameters are read flat, e.g.
`params.sigma_y`, or through the marker field, e.g. `params.cm_params.sigma_y`.

See also [`@dmg_params`](@ref), [`ConstitutiveParameters`](@ref), [`@params`](@ref).
"""
macro cm_params(model, params)
    macrocheck_input_model(model)
    macrocheck_input_model_params(params)
    return __model_params(:cm, model, params, __module__)
end

"""
    @dmg_params dmgmodel struct ModelParameters ... end

$(extension_api_note())

This is [`@cm_params`](@ref) for a damage model: the generated parameter type resolves the
`dmg_params::DamageParameters` marker field of a material. The standard fracture
parameters are declared this way. [`CriticalStretch`](@ref) owns `Gc` and `εc` through
the [`FractureParameters`](@ref) block:

```julia
Peridynamics.@dmg_params CriticalStretch struct CriticalStretchParameters
    @inherit FractureParameters
end
```

A custom damage model that wants the standard fracture keywords inherits the same block;
one that brings its own keywords declares them, and `material!` accepts them exactly when
the body's damage model reads them.

See also [`DamageParameters`](@ref), [`get_frac_params`](@ref), [`@params`](@ref).
"""
macro dmg_params(model, params)
    macrocheck_input_model(model)
    macrocheck_input_model_params(params)
    return __model_params(:dmg, model, params, __module__)
end

function __model_params(kind::Symbol, model, params_expr, mod::Module)
    local _name = params_expr.args[2]
    local _spec = get_param_decls(params_expr.args[3].args, mod)
    isempty(_spec.decls) && throw(ArgumentError(empty_model_params_msg(kind, _name)))
    check_model_param_block(kind, _spec)
    local _sim_float = params_use_sim_float(_spec)
    local _supertype = kind === :cm ? :(Peridynamics.AbstractConstitutiveParameters) :
                       :(Peridynamics.AbstractDamageParameters)
    local _struct = params_struct_expr(params_expr, _name, _supertype, _spec, _sim_float)
    local _convert = params_convert_expr(_name, _spec, _sim_float)
    local _fields_expr = quote
        function Peridynamics.param_fields_expr(::Base.Type{<:$(esc(_name))})
            return $(QuoteNode(_spec))
        end
    end
    local _logs = param_log_methods(_spec, _name)
    local _constructor = model_params_constructor_expr(model, _name, _spec, _sim_float,
                                                       mod)
    local _interface = model_params_interface_exprs(kind, model, _name, _spec, _sim_float)
    local _checks = quote
        Peridynamics.typecheck_model_params($(QuoteNode(kind)), $(esc(model)))
    end
    local _doc = quote
        Base.@__doc__ $(esc(_name))
    end
    return Expr(:block, _struct, _convert, _fields_expr, _logs..., _constructor,
                _interface..., _checks, _doc)
end

#=
The constructor mirrors the generated point parameter constructor: the model instance is
`model`, the material is `mat`, and every referenced name that is neither declared in the
block nor resolvable as a global is read from `mat_params`, the `NamedTuple` of the
material parameters declared above the marker field.
=#
function model_params_constructor_expr(model, name, spec, sim_float, mod::Module)
    body = param_constructor_body(spec, sim_float)
    args = [esc(d.name) for d in spec.decls]
    reads = material_level_reads(spec, mod)
    destructures = [Expr(:(=), esc(n),
                         :(Peridynamics.required_model_param($(esc(:mat_params)),
                                                             $(QuoteNode(n)),
                                                             $(esc(:model)))))
                    for n in reads]
    signature_args = (:($(esc(:model))::$(esc(model))),
                      :($(esc(:mat))::Peridynamics.AbstractMaterial),
                      :($(esc(:mat_params))::Base.NamedTuple),
                      :($(esc(:p))::Base.Dict{Base.Symbol,Base.Any}))
    sim_float || return quote
        function $(esc(name))($(signature_args...))
            $(destructures...)
            $(body...)
            return $(esc(name))($(args...))
        end
    end
    instantiation = Expr(:curly, esc(name), FLOAT_TYPE_PARAM)
    default_instantiation = Expr(:curly, esc(name), :(Peridynamics.default_float_type()))
    return quote
        function $(instantiation)($(signature_args...)) where {$(FLOAT_TYPE_PARAM)<:Base.Real}
            $(destructures...)
            $(body...)
            return $(instantiation)($(args...))
        end
        function $(esc(name))($(signature_args...))
            return $(default_instantiation)($(esc(:model)), $(esc(:mat)),
                                            $(esc(:mat_params)), $(esc(:p)))
        end
    end
end

function material_level_reads(spec::ParamFieldsSpec, mod::Module)
    declared = Set{Symbol}(d.name for d in spec.decls)
    reads = Symbol[]
    for decl in spec.decls
        for expr in (decl.default, decl.provider)
            isnothing(expr) && continue
            for n in referenced_names(expr)
                (n === :model || n === :mat || n === :p) && continue
                in(n, declared) && continue
                isnothing(resolve_in_mod_or_peridynamics(mod, n)) || continue
                in(n, reads) || push!(reads, n)
            end
        end
    end
    return reads
end

function model_params_interface_exprs(kind::Symbol, model, name, spec, sim_float)
    type_fn = kind === :cm ? :(Peridynamics.constitutive_param_type) :
              :(Peridynamics.damage_param_type)
    get_fn = kind === :cm ? :(Peridynamics.get_cm_params) : :(Peridynamics.get_dmg_params)
    kw_fn = kind === :cm ? :(Peridynamics.constitutive_param_kwargs) :
            :(Peridynamics.damage_param_kwargs)
    type_expr, get_expr, convert_expr = if sim_float
        (quote
             function $(type_fn)(::$(esc(model)),
                                 ::Base.Type{FT}=Peridynamics.default_float_type()) where {FT}
                 return $(Expr(:curly, esc(name), :FT))
             end
         end,
         quote
             function $(get_fn)(model::$(esc(model)), ::Base.Type{FT},
                                mat::Peridynamics.AbstractMaterial,
                                mat_params::Base.NamedTuple,
                                p::Base.Dict{Base.Symbol,Base.Any}) where {FT}
                 return $(Expr(:curly, esc(name), :FT))(model, mat, mat_params, p)
             end
         end,
         quote
             function Peridynamics.convert_nested_params(::Base.Type{FT_TO},
                                                         mp::$(esc(name))) where {FT_TO}
                 return $(Expr(:curly, esc(name), :FT_TO))(mp)
             end
         end)
    else
        (quote
             function $(type_fn)(::$(esc(model)),
                                 ::Base.Type=Peridynamics.default_float_type())
                 return $(esc(name))
             end
         end,
         quote
             function $(get_fn)(model::$(esc(model)), ::Base.Type,
                                mat::Peridynamics.AbstractMaterial,
                                mat_params::Base.NamedTuple,
                                p::Base.Dict{Base.Symbol,Base.Any})
                 return $(esc(name))(model, mat, mat_params, p)
             end
         end,
         quote
             Peridynamics.convert_nested_params(::Base.Type, mp::$(esc(name))) = mp
         end)
    end
    kw_expr = quote
        $(kw_fn)(::$(esc(model))) = $(Expr(:tuple, (QuoteNode(k) for k in spec.kwargs)...))
    end
    return (type_expr, get_expr, convert_expr, kw_expr)
end

function check_model_param_block(kind::Symbol, spec::ParamFieldsSpec)
    model = kind === :cm ? "constitutive model" : "damage model"
    for decl in spec.decls
        if is_model_param_decl(decl)
            msg = "a $(model) parameter block cannot carry the marker field "
            msg *= "`$(decl.name)`!\n"
            msg *= "  The markers belong to the point parameters of a material; a model "
            msg *= "cannot carry\n  the parameters of another model.\n"
            throw(ArgumentError(msg))
        end
    end
    return nothing
end

function typecheck_model_params(kind::Symbol, ::Type{Model}) where {Model}
    expected = kind === :cm ? AbstractConstitutiveModel : AbstractDamageModel
    Model <: expected && return nothing
    macro_name = kind === :cm ? "@cm_params" : "@dmg_params"
    msg = "`$(Model)` is not a subtype of `$(nameof(expected))`!\n"
    msg *= "  `$(macro_name)` declares the parameters of a $(nameof(expected)).\n"
    return throw(ArgumentError(msg))
end

function typecheck_model_params(kind::Symbol, model)
    return throw(ArgumentError("`$(model)` is not a type!\n"))
end

function empty_model_params_msg(kind::Symbol, name)
    macro_name = kind === :cm ? "@cm_params" : "@dmg_params"
    msg = "`$(macro_name) ... struct $(name)` declares no parameters!\n"
    msg *= "  A model without parameters simply does not use the macro.\n"
    return msg
end

function macrocheck_input_model(model)
    model isa Symbol && return nothing
    if model isa Expr && model.head === :.
        return nothing
    end
    msg = "argument `$(model)` is not a valid model type input!\n"
    return throw(ArgumentError(msg))
end

function macrocheck_input_model_params(params_expr)
    if !(params_expr isa Expr && params_expr.head === :struct)
        msg = "expected a `struct` definition, got `$(params_expr)`!\n"
        throw(ArgumentError(msg))
    end
    if !(params_expr.args[2] isa Symbol)
        msg = "the struct of a model parameter definition takes no supertype and no type "
        msg *= "parameters: `$(params_expr.args[2])`!\n"
        msg *= "  The macro derives both from the declarations.\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

