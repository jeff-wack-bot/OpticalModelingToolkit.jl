# Co-simulation components (Phase 3 of the RTSCosim plan): the minimal pieces
# a cavity length-control loop needs beyond the Phase 1 optics.
#
# Both follow the conventions of components.jl: plane-wave TEM₀₀ envelope as
# a length-2 real array on directed OpticalIn/OpticalOut sub-ports.

# -----------------------------------------------------------------------------
# Photodetector
# -----------------------------------------------------------------------------
#
# Non-loading power sensor: P = |a|² of the field entering `port_in`, exposed
# as an `[output = true]` signal variable. Terminates the beam (perfect
# absorber): its `port_out` emits zero, so it can be wired face-to-face like
# any other component.

@component function Photodetector(; name)
    systems = @named begin
        port_in = OpticalIn()
        port_out = OpticalOut()
    end

    # a guess, not a default: P is algebraically determined by the field, and
    # a default would overdetermine initialization
    @variables P(t) [output = true, description = "Detected power |a|²"]

    equations = Equation[
        P ~ port_in.a[1]^2 + port_in.a[2]^2,
        port_out.a[1] ~ 0,
        port_out.a[2] ~ 0,
    ]

    return System(equations, t, [P], []; name, systems, guesses = [P => 0.0])
end

# -----------------------------------------------------------------------------
# MovableMirror
# -----------------------------------------------------------------------------
#
# The `Mirror` scattering with a microscopic longitudinal position input `x`
# (an `[input = true]` signal, to be driven by a mechanical oscillator or a
# controller). Displacement by `x` toward the *front* side lengthens the
# back-side round trip and shortens the front-side one, so with the
# propagation convention `out = e^{ikL} in` of `Length`:
#
#     b_out.a = t·f_in.a + r·e^{+2ikx}·b_in.a
#     f_out.a = -r·e^{-2ikx}·f_in.a + t·b_in.a
#
# (transmission phase is position-independent to first order for a thin
# mirror). At x = 0 this reduces exactly to `Mirror`. The full nonlinear
# sin/cos dependence is kept — a length-control loop sweeps through fringes
# during acquisition.

@component function MovableMirror(; name, R = 0.99, T = nothing, λ = 1064e-9)
    Tval = T === nothing ? 1 - R : T
    pars = @parameters begin
        mr = sqrt(R), [description = "Amplitude reflectivity"]
        mt = sqrt(Tval), [description = "Amplitude transmissivity"]
        λ = λ, [description = "Carrier wavelength [m]"]
    end

    # a guess, not a default: x is driven by whatever the caller connects it
    # to, and a default would overdetermine initialization
    @variables x(t) [input = true,
        description = "Microscopic longitudinal position [m], positive toward the front side"]

    systems = @named begin
        f_in = OpticalIn()
        f_out = OpticalOut()
        b_in = OpticalIn()
        b_out = OpticalOut()
    end

    θ = 2 * (2π / λ) * x
    cθ = cos(θ)
    sθ = sin(θ)

    equations = Equation[
        # back-side reflection rotated by +θ
        b_out.a[1] ~ mt * f_in.a[1] + mr * (cθ * b_in.a[1] - sθ * b_in.a[2]),
        b_out.a[2] ~ mt * f_in.a[2] + mr * (sθ * b_in.a[1] + cθ * b_in.a[2]),
        # front-side reflection rotated by -θ
        f_out.a[1] ~ -mr * (cθ * f_in.a[1] + sθ * f_in.a[2]) + mt * b_in.a[1],
        f_out.a[2] ~ -mr * (-sθ * f_in.a[1] + cθ * f_in.a[2]) + mt * b_in.a[2],
    ]

    return System(equations, t, [x], pars; name, systems, guesses = [x => 0.0])
end
