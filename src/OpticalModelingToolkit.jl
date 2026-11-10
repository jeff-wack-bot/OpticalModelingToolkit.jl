module OpticalModelingToolkit

# OpticalModelingToolkit adds an optical physical domain to ModelingToolkit.

using LinearAlgebra: lu, issuccess, I
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SymbolicIndexingInterface
using SymbolicIndexingInterface: state_values, parameter_values
using Symbolics
using Symbolics: SymbolicUtils
import SciMLBase
using SimpleNonlinearSolve: SimpleNewtonRaphson

include("ports.jl")
export OpticalIn, OpticalOut, OpticalPort, link, expose

include("components.jl")
export Mirror, MovableMirror, Beamsplitter, FreeSpace
export Laser, LaserModulated, Vacuum, Photodetector
export PropagationRep, Algebraic, PadeDelay
export OpticalComponent, Scatterer, Propagator, Terminal

include("analysis.jl")
export linearize_proper, steady_state, frequency_response

include("scattering.jl")
export to_scattering

include("electronics.jl")
export ZPK

include("katscript.jl")
export parse_kat, build_kat

end # module OpticalModelingToolkit
