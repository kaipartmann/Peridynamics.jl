const IS_CI = haskey(ENV, "GITHUB_ACTIONS")
const LIVE_MODE = "LIVE_MODE" in ARGS
if LIVE_MODE
    using Revise
    Revise.revise()
end

using Peridynamics
using Documenter
using Literate

LIT_MD_OUT = joinpath(@__DIR__, "src", "generated")
# LIT_NB_OUT = joinpath(@__DIR__, "..", "notebooks") #TODO
rm(LIT_MD_OUT; recursive = true, force = true)
# rm(LIT_NB_OUT; recursive = true, force = true) #TODO

LIT_MANUAL_IN = [
    "howto_single_body_analysis.jl",
    "howto_contact_analysis.jl",
    "howto_pointclouds.jl",
    "howto_precracks.jl",
    "howto_matformulations.jl",
]
LIT_MANUAL_IN .= joinpath.(@__DIR__, "src", "literate", LIT_MANUAL_IN)
Literate.markdown.(LIT_MANUAL_IN, LIT_MD_OUT; credit=false)

LIT_TUTORIALS_IN = [
    "tutorial_tension_static.jl",
    "tutorial_tension_dynfrac.jl",
    "tutorial_tension_precrack.jl",
    "tutorial_wave_in_bar.jl",
]
LIT_TUTORIALS_IN .= joinpath.(@__DIR__, "src", "literate", LIT_TUTORIALS_IN)
Literate.markdown.(LIT_TUTORIALS_IN, LIT_MD_OUT; credit=false)
# Literate.notebook.(LIT_TUTORIALS_IN, LIT_NB_OUT; execute = IS_CI) #TODO

DocMeta.setdocmeta!(Peridynamics, :DocTestSetup, :(using Peridynamics); recursive=true)

makedocs(;
    modules = [Peridynamics],
    authors = "Kai Partmann",
    repo = "https://github.com/kaipartmann/Peridynamics.jl/blob/{commit}{path}#{line}",
    sitename = "Peridynamics.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://kaipartmann.github.io/Peridynamics.jl",
        edit_link = "main",
        assets = [joinpath("assets", "custom.css")],
        collapselevel = 1,
    ),
    draft = LIVE_MODE,
    pages = [
        "Home" => "index.md",
        "How-to guides" => [
            joinpath("generated", "howto_single_body_analysis.md"),
            joinpath("generated", "howto_contact_analysis.md"),
            joinpath("generated", "howto_pointclouds.md"),
            joinpath("generated", "howto_precracks.md"),
            joinpath("generated", "howto_matformulations.md"),
            "howto_visualization.md",
        ],
        "Tutorials" => [
            joinpath("generated", "tutorial_tension_static.md"),
            joinpath("generated", "tutorial_tension_dynfrac.md"),
            joinpath("generated", "tutorial_tension_precrack.md"),
            joinpath("generated", "tutorial_wave_in_bar.md"),
            "tutorial_logo.md",
        ],
        "API" => [
            "api_public.md",
            # "api_private.md", #TODO
        ]
    ],
)

if !LIVE_MODE
    deploydocs(;
        repo = "github.com/kaipartmann/Peridynamics.jl",
        devbranch = "main",
    )
end
