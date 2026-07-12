# Testing and Documentation System

## Structure

Every `.jl` file in `test/` (except `runtests.jl`) is a
[Literate.jl](https://github.com/fredrikekre/Literate.jl) source file. Each file
serves two purposes:

1. **As a test**: `runtests.jl` auto-discovers and includes every other `.jl`
   file in `test/`, wrapping each in a `@testset` named after the file.
2. **As a docs page**: `docs/make.jl` processes the same files with
   `Literate.markdown`, generating pages under `docs/src/generated/`. These
   appear in the "Tests" section of the docs.

Lines starting with `# ` become markdown exposition in the generated docs (use
`# ## Heading` for section headers and ```` ```math ```` blocks for display
equations). Code and `@test` assertions appear as executable code blocks.

## Prose vs. code comments

Literate classifies every line by its first non-whitespace character, ignoring
Julia's block structure:

- `# text` — becomes **markdown prose**. Only use at top level, *between*
  complete statements. A prose line in the middle of a `begin … end`,
  `function`, or `for` block will split the surrounding code block in the
  generated page and break execution.
- `## text` — stays a **code comment** inside the emitted code block. Use this
  for any explanatory comment that lives inside a `@testset`, function, or loop.

## Adding a new test file

Create a new `.jl` file in `test/`. It is picked up automatically by both the
test runner and the docs build — no changes to `runtests.jl` or `docs/make.jl`
are needed.

## Running tests

```
] test
```

The `test` target pulls in `ControlSystemsBase`, `NonlinearSolve`, `SciMLBase`,
`Symbolics`, `PythonCall`, and `Test` (see the `[targets]` section of the
top-level `Project.toml`). Finesse is installed by `CondaPkg.toml` via pip.

## Building docs

Full local build with live reload:

```julia
include("docs/serve.jl")
```

This calls `servedocs` with `skip_dir="docs/src/generated"` to avoid a rebuild
loop from Literate's generated output.

Targeted build of a single file to HTML (no Documenter required):

```julia
using Literate
Literate.html("test/test_propagation.jl", "output/")
```
