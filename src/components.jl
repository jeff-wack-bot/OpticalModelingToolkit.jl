# Phase 1a optical components.
#
# Plane-wave (TEM₀₀) optics with purely algebraic field relations. Every
# component is a pure constraint between its input and output sub-ports — no
# state, no time derivatives. The whole network compiles to a NonlinearSystem
# whose solution is the steady-state field amplitude at every node.
#
# Sign conventions follow plans/REFINED_PLAN.md §4 and the user-defined mirror
# scattering convention.

# -----------------------------------------------------------------------------
# Mirror
# -----------------------------------------------------------------------------
#
# A two-port lossless dielectric mirror with real-valued amplitude reflectivity
# `r` and amplitude transmissivity `t`, satisfying r² + t² = 1.
#
# The component has two physical sides ("front" / `f` and "back" / `b`), each
# split into directed sub-ports `*_in` (field entering) and `*_out` (field
# leaving). Scattering equations:
#
#     b_out.a = t * f_in.a + r * b_in.a
#     f_out.a = -r * f_in.a + t * b_in.a
#
# i.e. the scattering matrix is S = [[t, r], [-r, t]], which satisfies
# S†S = (t² + r²) I and is therefore unitary when t² + r² = 1.

@component function Mirror(; name, R = 0.99, T = nothing)
    Tval = T === nothing ? 1 - R : T
    # Parameters are named `mr`/`mt` (not `r`/`t`) to avoid shadowing the
    # MTK time variable `t` inside this function scope.
    pars = @parameters begin
        mr = sqrt(R), [description = "Amplitude reflectivity"]
        mt = sqrt(Tval), [description = "Amplitude transmissivity"]
    end

    systems = @named begin
        f_in = OpticalIn()
        f_out = OpticalOut()
        b_in = OpticalIn()
        b_out = OpticalOut()
    end

    equations = Equation[
        b_out.a[1] ~ mt * f_in.a[1] + mr * b_in.a[1],
        b_out.a[2] ~ mt * f_in.a[2] + mr * b_in.a[2],
        f_out.a[1] ~ -mr * f_in.a[1] + mt * b_in.a[1],
        f_out.a[2] ~ -mr * f_in.a[2] + mt * b_in.a[2],
    ]

    return System(equations, t, [], pars; name, systems)
end

# -----------------------------------------------------------------------------
# Length (algebraic phase shift)
# -----------------------------------------------------------------------------
#
# Phase 1a representation of free-space propagation: a purely algebraic phase
# rotation `out.a = exp(i k L) * in.a`. Length is unidirectional — to wire two
# mirrors face-to-face you instantiate one `Length` for the right-going field
# and a second one for the left-going field.
#
# In Phase 1b this same component will be promoted to a single-pole delay form
# `τ ȯ + o = exp(ikL) i` so that linearization yields a non-trivial state.

@component function Length(; name, L = 1.0, λ = 1064e-9)
    pars = @parameters begin
        L = L, [description = "Propagation length [m]"]
        λ = λ, [description = "Carrier wavelength [m]"]
    end

    systems = @named begin
        port_in = OpticalIn()
        port_out = OpticalOut()
    end

    # k * L = (2π/λ) * L, real-valued, time-independent.
    φ = 2π * L / λ
    cφ = cos(φ)
    sφ = sin(φ)

    equations = Equation[
        port_out.a[1] ~ cφ * port_in.a[1] - sφ * port_in.a[2],
        port_out.a[2] ~ sφ * port_in.a[1] + cφ * port_in.a[2],
    ]

    return System(equations, t, [], pars; name, systems)
end

# -----------------------------------------------------------------------------
# LengthDelay (single-pole approximation of free-space delay)
# -----------------------------------------------------------------------------
#
# Phase 1b representation of free-space propagation. The true delay
# `o(t) = exp(ikL) i(t - L/c)` is approximated by the first-order Padé form
#
#     τ ȯ + o = exp(ikL) i,    τ = L / c
#
# At DC (ω = 0) this matches the algebraic `Length` component exactly. At
# higher modulation frequencies the response rolls off as `1/(1 + iωτ)`,
# which is a one-pole approximation of the true `e^{-iωτ}` delay. Crucially,
# unlike the algebraic `Length`, this component introduces *state* — two real
# state variables per instance, one per quadrature — so the cavity has a
# non-trivial `(A, B, C, D)` after `linearize_symbolic`.

@component function LengthDelay(; name, L = 1.0, λ = 1064e-9, c = 299792458.0)
    pars = @parameters begin
        L = L, [description = "Propagation length [m]"]
        λ = λ, [description = "Carrier wavelength [m]"]
        c = c, [description = "Speed of light [m/s]"]
    end

    systems = @named begin
        port_in = OpticalIn()
        port_out = OpticalOut()
    end

    φ = 2π * L / λ
    cφ = cos(φ)
    sφ = sin(φ)
    τ = L / c

    equations = Equation[
        τ * D(port_out.a[1]) + port_out.a[1] ~ cφ * port_in.a[1] - sφ * port_in.a[2],
        τ * D(port_out.a[2]) + port_out.a[2] ~ sφ * port_in.a[1] + cφ * port_in.a[2],
    ]

    return System(equations, t, [], pars; name, systems)
end

# -----------------------------------------------------------------------------
# Laser (constant-amplitude source)
# -----------------------------------------------------------------------------
#
# Emits a fixed complex envelope `(A_re, A_im)` from its single output port and
# absorbs whatever comes back on its single input port. In Phase 1a the input
# port is a sink — the laser does not respond to back-action.

@component function Laser(; name, A = 1.0, ϕ = 0.0)
    pars = @parameters begin
        A = A, [description = "Output field amplitude"]
        ϕ = ϕ, [description = "Output field phase [rad]"]
    end

    systems = @named begin
        port_out = OpticalOut()
        port_in = OpticalIn()  # absorbs reflected light; no equation needed
    end

    equations = Equation[
        port_out.a[1] ~ A * cos(ϕ),
        port_out.a[2] ~ A * sin(ϕ),
    ]

    return System(equations, t, [], pars; name, systems)
end

# -----------------------------------------------------------------------------
# LaserModulated (constant carrier + small-signal amplitude modulation input)
# -----------------------------------------------------------------------------
#
# Same as `Laser`, but exposes a small-signal in-phase amplitude modulation
# `δA(t)` as a variable marked `[input = true]`. Defaulting to zero, this
# variable does nothing at the DC operating point — but it can be passed to
# `linearize_symbolic` as the input around which the system is linearized,
# letting us compute the cavity transfer function from laser amplitude to
# any field of interest.

@component function LaserModulated(; name, A = 1.0, ϕ = 0.0)
    pars = @parameters begin
        A = A, [description = "Carrier field amplitude"]
        ϕ = ϕ, [description = "Carrier field phase [rad]"]
    end

    @variables δA(t) = 0.0 [input = true,
        description = "Small-signal in-phase amplitude perturbation"]

    systems = @named begin
        port_out = OpticalOut()
        port_in = OpticalIn()  # absorbs reflected light
    end

    equations = Equation[
        port_out.a[1] ~ (A + δA) * cos(ϕ),
        port_out.a[2] ~ (A + δA) * sin(ϕ),
    ]

    return System(equations, t, [δA], pars; name, systems)
end

# -----------------------------------------------------------------------------
# Vacuum (zero-amplitude source / beam dump)
# -----------------------------------------------------------------------------
#
# Emits zero field on its output port. Used to terminate an unused side of a
# mirror so the network has no dangling directed connectors.

@component function Vacuum(; name)
    systems = @named begin
        port_out = OpticalOut()
        port_in = OpticalIn()
    end

    equations = Equation[
        port_out.a[1] ~ 0,
        port_out.a[2] ~ 0,
    ]

    return System(equations, t, [], []; name, systems)
end
