#=
The tables that say what a block exposes. They are rendered from the declarations that
`@params`, `@params_fields`, `@storage` and `@storage_fields` register, so the documentation
of a block cannot drift away from the block: the docstring of every block interpolates
`block_table`, the reference page is generated from it, and so is the REPL output of the
block itself.
=#

"""
    block_table(x)

$(extension_api_note())

Return the markdown table of what `x` exposes, where `x` is one of

| argument | table |
|:---|:---|
| a parameter block of [`@params_fields`](@ref) | its parameters and `material!` keywords |
| a point parameter type of [`@params`](@ref) | the same, for the whole type |
| a material | the table of its point parameters |
| a field block of [`@storage_fields`](@ref) | its storage fields |
| a storage type of [`@storage`](@ref) | the same, for the whole storage |

These are exactly the things [`@inherit`](@ref) accepts, so this answers "which names do I get
if I inherit this".

# Example
```julia-repl
julia> print(Peridynamics.block_table(Peridynamics.DiscretizationParameters))
| parameter | type | `material!` keyword | value | simulation log |
|:---|:---|:---|:---|:---|
| `δ` | simulation float | `horizon` | required | horizon |
| `rho` | simulation float | `rho` | required | density |
```
"""
function block_table end

function block_table(::Type{T}) where {T}
    spec = block_spec(T)
    isnothing(spec) && throw(ArgumentError(no_block_msg(T)))
    return block_table(spec)
end

block_table(mat::AbstractMaterial) = block_table(point_param_type(mat))

function no_block_msg(T)
    msg = "`$(T)` does not declare any parameters or storage fields!\n"
    msg *= "  Only what `@params`, `@params_fields`, `@storage` and `@storage_fields` "
    msg *= "define has a table, which is exactly what `@inherit` accepts.\n"
    return msg
end

#=
The registered declarations of anything that can be inherited from, or `nothing`. The
`param_fields_expr` and `storage_fields_expr` fallbacks throw for a type that declares
nothing, which is the right behavior for `@inherit` but not here, where the two are tried in
turn.
=#
function block_spec(::Type{T}) where {T}
    T <: AbstractPointParameterFields && return param_fields_expr(T)
    T <: AbstractPointParameters && return param_fields_expr(T)
    T <: AbstractConstitutiveParameters && return param_fields_expr(T)
    T <: AbstractDamageParameters && return param_fields_expr(T)
    T <: AbstractStorageFields && return storage_fields_expr(T)
    T <: AbstractStorage && return storage_fields_expr(T)
    # the nested states declare their fields the same way a storage does, so they render
    # the same table
    T <: AbstractConstitutiveState && return storage_fields_expr(T)
    T <: AbstractDamageState && return storage_fields_expr(T)
    return nothing
end

# --------------------------------------------------------------------------------------
# point parameters
# --------------------------------------------------------------------------------------

function block_table(spec::ParamFieldsSpec)
    msg = "| parameter | type | `material!` keyword | value | simulation log |\n"
    msg *= "|:---|:---|:---|:---|:---|\n"
    for decl in spec.decls
        msg *= "| `$(decl.name)` | $(param_type_msg(decl)) | $(param_kwarg_msg(decl)) | "
        msg *= "$(param_value_msg(decl)) | $(label_msg(decl.label)) |\n"
    end
    msg *= param_groups_msg(spec)
    isempty(spec.kwargs) && return msg
    keywords = join(("`$(k)`" for k in spec.kwargs), ", ")
    return msg * "\nKeywords of `material!`: $(keywords).\n"
end

#=
The groups are listed as they were declared, which is where the `material!` keywords and the
parameters they turn into meet: the keyword arguments of the printed call are the keywords,
the destructured names are the parameters.
=#
function param_groups_msg(spec::ParamFieldsSpec)
    groups = Vector{Pair{String,Vector{Symbol}}}()
    for decl in spec.decls
        is_provided(decl) || continue
        if !isempty(groups) && last(groups).first == decl.source
            push!(last(groups).second, decl.name)
        else
            push!(groups, decl.source => [decl.name])
        end
    end
    isempty(groups) && return ""
    msg = "\nComputed together:\n"
    for (source, names) in groups
        msg *= "- `(; $(join(names, ", "))) = $(source)`\n"
    end
    return msg
end

function param_type_msg(decl::ParamFieldDecl)
    decl.type === SimFloat && return "simulation float"
    return "`$(type_expr_string(decl.type))`"
end

param_kwarg_msg(decl::ParamFieldDecl) = decl.kwarg === :none ? "–" : "`$(decl.kwarg)`"

#=
The right-hand side is shown as it was written, which is what makes the table answer both
questions at once: which parameters a call supplies, and which `material!` keywords it reads,
because those are the keyword arguments of the very call that is printed.
=#
function param_value_msg(decl::ParamFieldDecl)
    is_cm_param_decl(decl) && return "owned by the constitutive model"
    is_dmg_param_decl(decl) && return "owned by the damage model"
    is_provided(decl) && return "from `$(first(split(decl.source, "(")))`"
    isnothing(decl.default) && return "required"
    return "`= $(decl.source)`"
end

# --------------------------------------------------------------------------------------
# storage fields
# --------------------------------------------------------------------------------------

function block_table(decls::Vector{StorageFieldDecl})
    msg = "| field | shape | entries | halo exchange |\n"
    msg *= "|:---|:---|:---|:---|\n"
    for decl in decls
        msg *= "| `$(decl.name)` | `$(type_msg(decl.type))` | $(entries_msg(decl)) | "
        msg *= "$(exchange_msg(decl)) |\n"
    end
    return msg
end

function entries_msg(decl::StorageFieldDecl)
    is_cm_state_decl(decl) && return "state of the constitutive model"
    is_dmg_state_decl(decl) && return "state of the damage model"
    decl.shape isa AbstractPointFieldShape && return "points"
    decl.shape isa AbstractBondFieldShape && return "bonds"
    decl.type === DofVector && return "degrees of freedom"
    return "–"
end

function exchange_msg(decl::StorageFieldDecl)
    decl.annotation === :lth && return "local → halo"
    decl.annotation === :htl && return "halo → local"
    return "–"
end

# --------------------------------------------------------------------------------------
# shared
# --------------------------------------------------------------------------------------

label_msg(label::AbstractString) = isempty(label) ? "–" : label

#=
Typing the name of a block at the REPL is the fastest way to ask what it exposes, so it
answers with its table instead of with its own name. The markdown is rendered, so the table
is aligned in a terminal and is a real table in the documentation. Only the blocks themselves
do this, not the point parameters and storages that `@params` and `@storage` generate: those
are printed as part of a `Body` and of error messages, where a table would be in the way.
=#
const BlockFields = Union{AbstractPointParameterFields,AbstractStorageFields}

function Base.show(io::IO, mime::MIME"text/plain", ::Type{T}) where {T<:BlockFields}
    return show_block(io, mime, T)
end

function Base.show(io::IO, mime::MIME"text/html", ::Type{T}) where {T<:BlockFields}
    return show_block(io, mime, T)
end

function show_block(io::IO, mime::MIME, ::Type{T}) where {T<:BlockFields}
    spec = try
        block_spec(T)
    catch
        nothing
    end
    # `print` rather than `show(io, mime, T)`, which would be this method again
    isnothing(spec) && return print(io, T)
    msg = "**`$(T)`** – a block you can `@inherit`:\n\n" * block_table(spec)
    show(io, mime, Markdown.parse(msg))
    return nothing
end

# the tables are read by people who write `PointVector{Float64}`, not `Peridynamics.PointVector{Float64}`
type_msg(type) = replace(string(type), "Peridynamics." => "")
