using MTKOptical
using Documenter

DocMeta.setdocmeta!(MTKOptical, :DocTestSetup, :(using MTKOptical); recursive=true)

makedocs(;
    modules=[MTKOptical],
    authors="Jeff Wack <jeffwack111@gmail.com> and contributors",
    sitename="MTKOptical.jl",
    format=Documenter.HTML(;
        canonical="https://jeff.github.io/MTKOptical.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jeff/MTKOptical.jl",
    devbranch="main",
)
