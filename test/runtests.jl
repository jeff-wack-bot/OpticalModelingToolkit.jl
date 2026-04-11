using OpticalModelingToolkit
using Test

@testset "OpticalModelingToolkit.jl" begin

    @testset "Finesse Integration Tests" begin
        include("test_finesse_simple.jl")
    end

    @testset "Phase 1a Fabry-Perot" begin
        include("test_phase1a_fabry_perot.jl")
    end

    @testset "Phase 1b Fabry-Perot AC" begin
        include("test_phase1b_fabry_perot_ac.jl")
    end

end
