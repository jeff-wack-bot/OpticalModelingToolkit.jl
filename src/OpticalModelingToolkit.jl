module OpticalModelingToolkit

# OpticalModelingToolkit adds an optical physical domain to ModelingToolkit.
# See plans/REFINED_PLAN.md for the full design.

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary
using Symbolics

# PythonCall is used by the Finesse exporter (Phase 2+) and by the validation
# tests, which build an equivalent Finesse model in Python and compare results
# against the MTK model. Loading PythonCall here triggers CondaPkg to
# materialize the project-local Python env on first use.
using PythonCall

include("connectors.jl")
export OpticalIn, OpticalOut

include("components.jl")
export Mirror, Length, LengthDelay, Laser, LaserModulated, Vacuum

end # module OpticalModelingToolkit
