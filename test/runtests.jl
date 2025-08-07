using MTKOptical
using Test

@testset "MTKOptical.jl" begin
    
    @testset "Basic Python-Julia Interoperability" begin
        include("test_basic_interop.jl")
    end
    
    @testset "Finesse Integration Tests" begin
        include("test_finesse_simple.jl")
    end
    
end
