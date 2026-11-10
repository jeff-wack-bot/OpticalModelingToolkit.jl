module OpticalModelingToolkit

# OpticalModelingToolkit adds an optical physical domain to ModelingToolkit.
# See plans/REFINED_PLAN.md for the full design.

using LinearAlgebra: lu, issuccess
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SymbolicIndexingInterface: state_values, parameter_values
using Symbolics
using Symbolics: SymbolicUtils


include("connectors.jl")
export OpticalIn, OpticalOut

include("components.jl")
export Mirror, FreeSpace, Laser, LaserModulated, Vacuum
export PropagationRep, Algebraic, PadeDelay, OpticalComponent

include("cosim_components.jl")
export Photodetector, MovableMirror

include("analysis.jl")
export linearize_proper

include("scattering.jl")
export to_scattering

end # module OpticalModelingToolkit
