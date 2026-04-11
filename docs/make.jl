using OpticalModelingToolkit
using Documenter

DocMeta.setdocmeta!(OpticalModelingToolkit, :DocTestSetup, :(using OpticalModelingToolkit); recursive=true)

makedocs(;
    modules=[OpticalModelingToolkit],
    authors="Jeff Wack <jeffwack111@gmail.com> and contributors",
    sitename="OpticalModelingToolkit.jl",
    format=Documenter.HTML(;
        canonical="https://jeff.github.io/OpticalModelingToolkit.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Project Plan" => "PROJECT_PLAN.md",
        "Python-Julia Interoperability" => "PYTHON_JULIA_INTEROP.md",
        "Finesse Analysis" => "FINESSE_ANALYSIS.md",
    ],
)

deploydocs(;
    repo="github.com/jeff/OpticalModelingToolkit.jl",
    devbranch="main",
)
