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
LIT_NB_OUT = joinpath(@__DIR__, "..", "notebooks")
rm(LIT_MD_OUT; recursive = true, force = true)
rm(LIT_NB_OUT; recursive = true, force = true)

LIT_MANUAL_IN = ["howto_spatial_discretization.jl",
                 "howto_matmodels.jl",
                 "howto_temporal_discretization.jl",
                 "howto_conditions.jl",
                 "howto_jobs.jl",
                 "howto_workflow.jl",
                 "howto_vtk_reader.jl"]
LIT_MANUAL_IN .= joinpath.(@__DIR__, "src", "literate", LIT_MANUAL_IN)
Literate.markdown.(LIT_MANUAL_IN, LIT_MD_OUT)

LIT_TUTORIALS_IN = []
LIT_TUTORIALS_IN .= joinpath.(@__DIR__, "src", "literate", LIT_TUTORIALS_IN)
Literate.markdown.(LIT_TUTORIALS_IN, LIT_MD_OUT)
Literate.notebook.(LIT_TUTORIALS_IN, LIT_NB_OUT; execute = IS_CI)

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
        assets = String[],
        collapselevel = 1,
    ),
    draft = LIVE_MODE,
    pages = [
        "Home" => "index.md",
        "quickstart.md",
        "theory.md",
        "How-to guides" => [
            joinpath("generated", "howto_spatial_discretization.md"),
            joinpath("generated", "howto_matmodels.md"),
            joinpath("generated", "howto_temporal_discretization.md"),
            joinpath("generated", "howto_conditions.md"),
            joinpath("generated", "howto_jobs.md"),
            joinpath("generated", "howto_workflow.md"),
            "howto_visualization.md",
            joinpath("generated", "howto_vtk_reader.md"),
        ],
        "Tutorials" => [
            "tutorial_tensiletest.md",
            "tutorial_crackedplateundertension.md",
            "tutorial_logo.md",
        ],
        "API" => [
            "api_public.md",
            "api_private.md",
        ]
    ],
)

if !LIVE_MODE
    deploydocs(;
        repo = "github.com/kaipartmann/Peridynamics.jl",
        devbranch = "main",
    )
end
