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

LIT_MANUAL_IN = ["spatial_discretization.jl",
                 "matmodels.jl",
                 "temporal_discretization.jl",
                 "conditions.jl",
                 "jobs.jl",
                 "workflow.jl",
                 "vtk_reader.jl"]
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
        "Manual" => [
            joinpath("generated", "spatial_discretization.md"),
            joinpath("generated", "matmodels.md"),
            joinpath("generated", "temporal_discretization.md"),
            joinpath("generated", "conditions.md"),
            joinpath("generated", "jobs.md"),
            joinpath("generated", "workflow.md"),
            "visualization.md",
            joinpath("generated", "vtk_reader.md"),
        ],
        "Examples" => [
            "tensiletest.md",
            "crackedplateundertension.md",
            "logo.md",
        ],
        "API" => [
            "public_api.md",
            "private_api.md",
        ]
    ],
)

if !LIVE_MODE
    deploydocs(;
        repo = "github.com/kaipartmann/Peridynamics.jl",
        devbranch = "main",
    )
end
