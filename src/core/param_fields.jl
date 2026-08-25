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
- `type::Any`: Declared type of the parameter, or [`SimFloat`](@ref) if the parameter
    follows the float type of the simulation.
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

@inline is_provided(d::ParamFieldDecl) = !isnothing(d.provider)
@inline is_derived(d::ParamFieldDecl) = isnothing(d.provider) && d.kwarg === :none

@inline function param_field_type(d::ParamFieldDecl)
    return d.type === SimFloat ? FLOAT_TYPE_PARAM : d.type
end

@inline function params_use_sim_float(spec::ParamFieldsSpec)
    return any(d -> d.type === SimFloat, spec.decls)
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
@derived (; Gc, εc) = get_frac_params(mat.dmgmodel, δ, K; Gc, epsilon_c)
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
logged; the second declares and labels in one line. Without this annotation the parameter
does not appear in the log, which is the default of `log_param_property`.

```julia
@log "yield stress" sigma_y = Inf
@log "bond constant" @derived bc = 18 * K / (π * δ^4)

@derived (; Gc, εc) = get_frac_params(mat.dmgmodel, δ, K; Gc, epsilon_c)
@log "critical stretch" εc
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
[`ElasticParameters`](@ref) and [`FractureParameters`](@ref) (combined in
[`StandardParameters`](@ref)), [`BondHorizonParameters`](@ref) and
[`InteractionParameters`](@ref), each documented with the parameters it exposes.
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
            push!(decls, ParamFieldDecl(member, SimFloat, :none, nothing, provider, source,
                                        ""))
        elseif member isa Expr && member.head === :(::) && member.args[1] isa Symbol
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
    msg *= "        @derived (; Gc, εc) = get_frac_params(mat.dmgmodel, δ, K; Gc, "
    msg *= "epsilon_c)\n"
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
    msg *= "        get_frac_params(mat.dmgmodel, δ, K; Gc, epsilon_c)\n"
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
    msg *= "        @derived (; Gc, εc) = get_frac_params(mat.dmgmodel, δ, K; Gc, "
    msg *= "epsilon_c)\n"
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

function parse_param_field(expr, mod::Module)
    default = nothing
    decl = expr
    if decl isa Expr && decl.head === :(=)
        default = decl.args[2]
        decl = decl.args[1]
    end
    if decl isa Symbol
        # an untyped parameter follows the float type of the simulation
        return (; name=decl, type=SimFloat, default)
    end
    if decl isa Expr && decl.head === :(::) && decl.args[1] isa Symbol
        return (; name=decl.args[1], type=resolve_param_type(mod, decl.args[2]), default)
    end
    return throw(ArgumentError("unexpected parameter declaration: $expr\n"))
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
    resolved = resolve_in_mod_or_peridynamics(mod, type_expr)
    isa(resolved, Type) && return resolved
    msg = "cannot resolve the type `$(type_expr)` of a point parameter!\n"
    msg *= "  It has to be resolvable in `$(mod)` or in `Peridynamics`. Declare the "
    msg *= "parameter without a type to let it follow the float type of the simulation.\n"
    return throw(ArgumentError(msg))
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
    d.type === SimFloat || (msg *= "::$(d.type)")
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
function param_constructor_body(spec::ParamFieldsSpec)
    body = Any[]
    provider, provider_var = nothing, nothing
    for decl in spec.decls
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
        msg *= log_param_property(Val(key), param; indentation)
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
