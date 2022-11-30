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
        "Manual" => [
            "manual/manual.md",
        ],
        "Tutorials" => [
            "tutorials/tensiletest.md",
            "tutorials/crackedplateundertension.md",
            "tutorials/logo.md",
            "tutorials/paraviewtutorial.md",
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
