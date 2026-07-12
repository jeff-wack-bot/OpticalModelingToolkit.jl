# Loaded when Etalon is present: wrap `to_scattering` output in Etalon's
# ScatteringSystem, inheriting its exact transcendental `freqresp` and
# `pade(sys; order)`. The interface is one-directional (OMTK → struct);
# neither package touches the other's internals.

module OpticalModelingToolkitEtalonExt

using OpticalModelingToolkit: to_scattering
using ModelingToolkit: System
using Etalon: Etalon

"""
    Etalon.ScatteringSystem(model::System)

Lower a composed OpticalModelingToolkit model to Etalon's exact transcendental
scattering representation — see [`to_scattering`](@ref).
"""
Etalon.ScatteringSystem(model::System) = Etalon.ScatteringSystem(to_scattering(model)...)

end # module
