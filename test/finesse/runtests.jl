# Cross-validation of OpticalModelingToolkit against Finesse.jl, a pure-Julia
# port of FINESSE 3 (https://git.mccullerlab.com/software/finesse.jl).
#
# This is a separate environment, not part of `Pkg.test()`: Finesse.jl and its
# dependency SignalFlowLU are unregistered, and OMT should not depend on them.
# Run from the repository root:
#
#     julia --project=test/finesse -e 'using Pkg; Pkg.instantiate()'
#     julia --project=test/finesse test/finesse/runtests.jl
#
# The OMT suite (`test/runtests.jl`) checks the package against analytic
# references and against itself. This one checks it against an independent
# frequency-domain simulator that was written from a different codebase, on the
# same networks: DC fields at every port, exact sideband transfer functions,
# the linearized laser-to-photodetector response, and the optical spring.

using Test
using LinearAlgebra
using OpticalModelingToolkit
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using NonlinearSolve
using SciMLBase
using ControlSystemsBase: ss, dcgain, freqresp

include("common.jl")

@testset "OpticalModelingToolkit vs Finesse.jl" begin
    for f in ("test_conventions.jl", "test_dc_fields.jl", "test_sidebands.jl", "test_signals.jl", "test_kat_models.jl")
        @testset "$f" begin
            include(f)
        end
    end
end
