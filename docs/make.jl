using OpticalModelingToolkit
using Documenter
using Literate

DocMeta.setdocmeta!(OpticalModelingToolkit, :DocTestSetup, :(using OpticalModelingToolkit); recursive = true)

# Every `.jl` file in `test/` (except `runtests.jl`) is a Literate.jl source
# file. Render each one as an Examples page, titled by its own heading.
tests_dir = joinpath(@__DIR__, "..", "test")
generated_dir = joinpath(@__DIR__, "src", "generated")
for testfile in filter(f -> endswith(f, ".jl") && f != "runtests.jl", readdir(tests_dir))
    Literate.markdown(joinpath(tests_dir, testfile), generated_dir; flavor = Literate.DocumenterFlavor())
end
examples = [joinpath("generated", f) for f in readdir(generated_dir) if endswith(f, ".md")]

makedocs(;
    modules = [OpticalModelingToolkit],
    authors = "Jeff Wack <jeffwack111@gmail.com> and contributors",
    sitename = "OpticalModelingToolkit.jl",
    format = Documenter.HTML(;
        canonical = "https://jeff-wack-bot.github.io/OpticalModelingToolkit.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Design" => "design.md",
        "Solver strategy" => "solvers.md",
        "Examples" => examples,
        "API" => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/jeff-wack-bot/OpticalModelingToolkit.jl",
    devbranch = "main",
    push_preview = true,
)
