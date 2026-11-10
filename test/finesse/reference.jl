# Write `models/reference/*.txt`: Finesse.jl's responses of the shared models.
#
#     julia --project=test/finesse test/finesse/reference.jl
#
# Run it after changing a model; `test_kat_models.jl` fails until you do.

using Test
using LinearAlgebra
using OpticalModelingToolkit
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using NonlinearSolve
using SciMLBase

include("common.jl")
include("kat_models.jl")

mkpath(joinpath(MODELS_DIR, "reference"))
for (model, pairs) in REFERENCES
    write_reference(model, pairs, finesse_references(model, pairs))
    println("wrote ", reference_path(model))
end
