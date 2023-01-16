using Peridynamics
using Documenter

DocMeta.setdocmeta!(Peridynamics, :DocTestSetup, :(using Peridynamics); recursive=true)

makedocs(;
    modules=[Peridynamics],
    authors="Kai Friebertshäuser",
    repo="https://github.com/kfrb/Peridynamics.jl/blob/{commit}{path}#{line}",
    sitename="Peridynamics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kfrb.github.io/Peridynamics.jl",
        edit_link="main",
        assets=String[],
        collapselevel=1,
    ),
    pages=[
        "Home" => "index.md",
        "quickstart.md",
        "Manual" => [
            "manual/spatial_discretization.md",
            "manual/matmodels.md",
            "manual/temporal_discretization.md",
            "manual/conditions.md",
            "manual/jobs.md",
            "manual/workflow.md",
            "manual/visualization.md",
            "manual/vtk_reader.md",
        ],
        "Examples" => [
            "examples/tensiletest.md",
            "examples/crackedplateundertension.md",
            "examples/logo.md",
        ],
        "API" => [
            "api/public_api.md",
            "api/private_api.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/kfrb/Peridynamics.jl",
    devbranch="main",
)
