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
        "manual.md",
        "Tutorials" => [
            "tensiletest.md",
            "crackedplateundertension.md",
            "logo.md",
            "paraviewtutorial.md",
        ],
        "API" => [
            "public_api.md",
            "private_api.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/kfrb/Peridynamics.jl",
    devbranch="main",
)
