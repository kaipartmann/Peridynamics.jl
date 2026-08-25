const IS_CI = haskey(ENV, "GITHUB_ACTIONS")
const LIVE_MODE = "LIVE_MODE" in ARGS

# The API reference is split by tier: `extension_api_reference.md` lists the names declared
# `public` in `src/public_api.jl`, and `internals.md` picks up everything else with
# `@autodocs Public = false`. Both depend on the `public` keyword, which exists from
# Julia 1.11 on. On 1.10 every public name would land on both pages and Documenter would
# fail on a duplicated docstring.
if VERSION < v"1.11"
    error("""
          building the documentation needs Julia 1.11 or newer, found $(VERSION).

          The API reference is split with `checkdocs = :public` and `@autodocs Public = ...`,
          which both need the `public` keyword. The package itself still supports 1.10.
          """)
end
if LIVE_MODE
    using Revise
    Revise.revise()
end

using Peridynamics
using Documenter
using Literate
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"), style=:alpha)

LIT_MD_OUT = joinpath(@__DIR__, "src", "generated")
# LIT_NB_OUT = joinpath(@__DIR__, "..", "notebooks") #TODO
rm(LIT_MD_OUT; recursive = true, force = true)
# rm(LIT_NB_OUT; recursive = true, force = true) #TODO

# LIT_MANUAL_IN = [
    # "howto_single_body_analysis.jl",
    # "howto_contact_analysis.jl",
    # "howto_pointclouds.jl",
    # "howto_precracks.jl",
    # "howto_matformulations.jl",
# ]
# LIT_MANUAL_IN .= joinpath.(@__DIR__, "src", "literate", LIT_MANUAL_IN)
# Literate.markdown.(LIT_MANUAL_IN, LIT_MD_OUT; credit=false)

LIT_TUTORIALS_IN = [
    "tutorial_tension_static.jl",
    "tutorial_tension_dynfrac.jl",
    "tutorial_tension_precrack.jl",
    "tutorial_wave_in_bar.jl",
    "tutorial_kalthoff-winkler_dynfrac.jl",
    "tutorial_logo.jl",
    "tutorial_cylinder.jl",
    "tutorial_wave_interface.jl",
    "tutorial_brazilian_test.jl",
    "tutorial_custom_material.jl",
    "tutorial_custom_damage_model.jl",
    "tutorial_custom_constitutive_model.jl",
]
LIT_TUTORIALS_IN .= joinpath.(@__DIR__, "src", "literate", LIT_TUTORIALS_IN)
Literate.markdown.(LIT_TUTORIALS_IN, LIT_MD_OUT; credit=false)
# Literate.notebook.(LIT_TUTORIALS_IN, LIT_NB_OUT; execute = IS_CI) #TODO

#=
"Blocks you can inherit" is generated rather than written, so that a block appears on it by
existing. It is the index of the "Blocks" section of the extension API reference, where the
docstrings live, and it needs nothing but the lists of blocks, which are the public subtypes
of the marker supertypes.
=#
using InteractiveUtils: subtypes

function public_blocks(supertype)
    blocks = filter(subtypes(supertype)) do T
        return parentmodule(T) === Peridynamics && Base.ispublic(Peridynamics, nameof(T))
    end
    return sort!(blocks; by=nameof)
end

function block_links(supertype)
    return join(("- [`$(nameof(T))`](@ref Peridynamics.$(nameof(T)))"
                 for T in public_blocks(supertype)), "\n")
end

function write_inheritable_blocks_page(path)
    msg = """
    # Blocks you can inherit

    [`@inherit`](@ref Peridynamics.@inherit) includes all declarations of another block into
    a [`@params`](@ref Peridynamics.@params) or [`@storage`](@ref Peridynamics.@storage)
    definition. This page lists everything the package ships that can be inherited. What a
    block exposes is on its reference entry, and the same table is printed by
    `Peridynamics.block_table(Block)` and by typing the name of a block at the REPL.

    Two `@inherit`s may contribute the same name only if they declare it identically, and a
    declaration in the body overrides an inherited one in place.

    ## Point parameter blocks

    The parameters a `@derived` right-hand side can read once the block is inherited, and
    the keywords [`material!`](@ref) then accepts.

    $(block_links(Peridynamics.AbstractPointParameterFields))

    ## Point parameters of the shipped materials

    A point parameter type can be inherited like a block, and it can also be used as it is
    with the second form of `@params`, e.g. `Peridynamics.@params MyMaterial BBPointParameters`.

    $(block_links(Peridynamics.AbstractPointParameters))

    ## Storage field blocks

    A storage needs the block of the time solver it is used with. The others follow from the
    system and the material family.

    $(block_links(Peridynamics.AbstractStorageFields))

    ## Storages of the shipped materials

    A material that builds on a family of this package inherits the storage of that family
    instead of listing its fields again, e.g. `@inherit RKCStorage`.

    $(block_links(Peridynamics.AbstractStorage))

    ## Blocks of your own

    ```julia
    Peridynamics.@params_fields PlasticityParameters begin
        @log "initial yield stress" @kwarg sigma_y σy = Inf
        @log "hardening modulus" @kwarg hardening Hiso = 0.0
    end
    ```

    Interpolate `\$(Peridynamics.block_table(PlasticityParameters))` into the docstring of a
    block of your own to give it the same table.
    """
    mkpath(dirname(path))
    write(path, msg)
    return nothing
end

write_inheritable_blocks_page(joinpath(LIT_MD_OUT, "inheritable_blocks.md"))

DocMeta.setdocmeta!(Peridynamics, :DocTestSetup, :(using Peridynamics); recursive=true)

makedocs(;
    plugins = [bib],
    modules = [Peridynamics],
    authors = "Kai Partmann",
    # a `Remote` rather than a URL template, so that Documenter can also build the navbar
    # link to the repository root
    repo = Remotes.GitHub("kaipartmann", "Peridynamics.jl"),
    sitename = "Peridynamics.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://kaipartmann.github.io/Peridynamics.jl",
        edit_link = "main",
        assets = [joinpath("assets", "custom.css")],
        collapselevel = 1,
        # The API reference pages are lists of docstrings and are inevitably large, and
        # `internals.md` is one `@autodocs` block over the whole package. Exempt exactly
        # these rather than raising the threshold for the prose pages, where a size warning
        # still catches something worth knowing.
        size_threshold_ignore = ["internals.md", "public_api_reference.md",
                                 "extension_api_reference.md"],
        # the search index covers every internal docstring through `internals.md`, so it
        # is larger than the default warn limit of 500 KiB; the hard limit stays at 1 MiB
        search_size_threshold_warn = 768 * 2^10,
    ),
    draft = LIVE_MODE,
    pages = [
        "Home" => "index.md",
        "Explanations" => [
            "expl_general_pd.md",
            "expl_bondbased.md",
            "expl_osbased.md",
            "expl_nosbased.md",
            "expl_continuumbased.md",
            "expl_damage.md",
            "expl_references.md",
        ],
        "How-to guides" => [
            "howto_mpi.md",
            "howto_visualization.md",
            "howto_parameter-study.md",
        ],
        "Tutorials" => [
            joinpath("generated", "tutorial_tension_static.md"),
            joinpath("generated", "tutorial_tension_dynfrac.md"),
            joinpath("generated", "tutorial_tension_precrack.md"),
            joinpath("generated", "tutorial_wave_in_bar.md"),
            joinpath("generated", "tutorial_wave_interface.md"),
            joinpath("generated", "tutorial_kalthoff-winkler_dynfrac.md"),
            joinpath("generated", "tutorial_logo.md"),
            joinpath("generated", "tutorial_cylinder.md"),
            joinpath("generated", "tutorial_brazilian_test.md"),
        ],
        "Development" => [
            joinpath("generated", "tutorial_custom_material.md"),
            joinpath("generated", "tutorial_custom_damage_model.md"),
            joinpath("generated", "tutorial_custom_constitutive_model.md"),
            "dev_systems.md",
            "dev_materials.md",
            joinpath("generated", "inheritable_blocks.md"),
            "dev_solvers.md",
            "dev_multithreading_mpi.md",
        ],
        "API Reference" => [
            "api_stability.md",
            "public_api_reference.md",
            "extension_api_reference.md",
            "internals.md",
        ]
    ],

    # Every exported and every `public` name must appear on one of the curated reference
    # pages. Internal docstrings are swept up by the `@autodocs` block in `internals.md`
    # instead of having to be listed by hand, which is what used to make this build drift.
    checkdocs = :public,
)

if !LIVE_MODE
    deploydocs(;
        repo = "github.com/kaipartmann/Peridynamics.jl",
        devbranch = "main",
    )
end
