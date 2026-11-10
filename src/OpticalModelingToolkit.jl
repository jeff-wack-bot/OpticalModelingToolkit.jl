module OpticalModelingToolkit

# OpticalModelingToolkit adds an optical physical domain to ModelingToolkit.

using LinearAlgebra: lu, issuccess
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SymbolicIndexingInterface: state_values, parameter_values
using Symbolics
using Symbolics: SymbolicUtils

include("ports.jl")
export OpticalIn, OpticalOut, OpticalPort, link

include("components.jl")
export Mirror, MovableMirror, Beamsplitter, FreeSpace
export Laser, LaserModulated, Vacuum, Photodetector
export PropagationRep, Algebraic, PadeDelay
export OpticalComponent, Scatterer, Propagator, Terminal

include("analysis.jl")
export linearize_proper

include("scattering.jl")
export to_scattering

end # module OpticalModelingToolkit
