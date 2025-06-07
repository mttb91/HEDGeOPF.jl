using HEDGeOPF
using Documenter

makedocs(;
    modules=[HEDGeOPF],
    authors="Matteo Baù",
    sitename="HEDGeOPF",
    format=Documenter.HTML(;
        canonical="https://mttb91.github.io/HEDGeOPF.jl",
        edit_link="main",
        mathengine = Documenter.MathJax(),
        assets=String[],
    ),
    checkdocs = :none,
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "quickstartguide.md",
            "Simulation Configuration" => "configuration.md",
            "Dataset Format" => "dataset.md"
        ]
    ],
)

deploydocs(;
    repo="github.com/mttb91/HEDGeOPF.jl",
    devbranch="main",
)
