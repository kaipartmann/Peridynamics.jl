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
| a parameter block of `@params_fields` | its parameters and `material!` keywords |
| a point parameter type of `@params` | the same, for the whole type |
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
    T <: AbstractStorageFields && return storage_fields_expr(T)
    T <: AbstractStorage && return storage_fields_expr(T)
    return nothing
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

#=
Typing the name of a block at the REPL is the fastest way to ask what it exposes, so it
answers with its table instead of with its own name. The markdown is rendered, so the table
is aligned in a terminal and is a real table in the documentation. Only the blocks themselves
do this, not the point parameters and storages that `@params` and `@storage` generate: those
are printed as part of a `Body` and of error messages, where a table would be in the way.
=#
const BlockFields = Union{AbstractStorageFields}

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
