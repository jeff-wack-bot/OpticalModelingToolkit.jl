using OpticalModelingToolkit
using Test

# Every `.jl` file in this directory (except this one) is a Literate.jl source
# file that doubles as a test. Auto-discover them and wrap each in a `@testset`
# named after the file. Adding a new test file requires no changes here.
@testset "OpticalModelingToolkit.jl" begin
    for f in filter(f -> endswith(f, ".jl") && f != "runtests.jl", readdir(@__DIR__))
        @testset "$f" begin
            include(f)
        end
    end
end
