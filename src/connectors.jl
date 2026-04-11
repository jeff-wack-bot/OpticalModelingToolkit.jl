# Optical connectors.
#
# Each Finesse-style optical port is split into two directed sub-ports: one
# carrying the field leaving the component (`OpticalOut`) and one carrying the
# field entering it (`OpticalIn`). Two components are wired face-to-face with
# two `connect` calls — one for each direction. See plans/REFINED_PLAN.md §4.
#
# The carried quantity is the slowly-varying complex envelope amplitude `a` of
# the optical field at the carrier frequency. MTK Standard Library is real-
# valued, so we represent the complex amplitude as a length-2 array
# `[real, imag]`. Modeled directly on `RealInputArray`/`RealOutputArray` from
# `ModelingToolkitStandardLibrary.Blocks` (utils.jl:35-119).

@connector function OpticalIn(; name, guess = zeros(2))
    @variables a(t)[1:2] [
        input = true,
        description = "Complex envelope amplitude (real, imag) entering the component"
    ]
    System(Equation[], t, [a], []; name, guesses = [a => guess])
end

@connector function OpticalOut(; name, guess = zeros(2))
    @variables a(t)[1:2] [
        output = true,
        description = "Complex envelope amplitude (real, imag) leaving the component"
    ]
    System(Equation[], t, [a], []; name, guesses = [a => guess])
end
