using OpticalModelingToolkit
using Documenter
using Literate

DocMeta.setdocmeta!(OpticalModelingToolkit, :DocTestSetup, :(using OpticalModelingToolkit); recursive=true)

# Every `.jl` file in `test/` (except `runtests.jl`) is a Literate.jl source
# file. Generate a docs page from each one.
tests_dir = joinpath(@__DIR__, "..", "test")
generated_dir = joinpath(@__DIR__, "src", "generated")
for testfile in filter(f -> endswith(f, ".jl") && f != "runtests.jl", readdir(tests_dir))
    Literate.markdown(joinpath(tests_dir, testfile), generated_dir; flavor=Literate.DocumenterFlavor())
end
test_pages = [splitext(f)[1] => joinpath("generated", f)
                 for f in filter(f -> endswith(f, ".md"), readdir(generated_dir))]

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
        "The Fabry-Perot cavity" => "fabry_perot.md",
        "Solver strategy" => "solvers.md",
        "Tests" => test_pages,
        "Testing Guide" => "TESTING.md",
    ],
)

deploydocs(;
    repo="github.com/jeff/OpticalModelingToolkit.jl",
    devbranch="main",
)
