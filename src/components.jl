# Optical components.
#
# Plane-wave (TEM₀₀) optics. Static components (Mirror, Laser, Vacuum) are pure
# constraints between their input and output sub-ports; free-space propagation
# is the one place dynamics can enter, and *how* it is represented is an
# analysis-time choice threaded in through `FreeSpace(rep = ...)` — never an
# element identity (plans/PLAN_rationalize_late.md §2).
#
# Sign conventions follow plans/REFINED_PLAN.md §4 and the user-defined mirror
# scattering convention.

# -----------------------------------------------------------------------------
# Component metadata
# -----------------------------------------------------------------------------
#
# Every component stamps its `System` metadata with `OpticalComponent =>
# (kind = ..., <numeric constructor parameters>)`. `to_scattering`
# (scattering.jl) walks this metadata to lower a composed model to the
# canonical scattering form without pattern-matching on equations.

"""
    OpticalComponent

`System`-metadata key under which every OpticalModelingToolkit component
records its kind and numeric constructor parameters, e.g.
`(kind = :Mirror, r = 0.995, t = 0.1)`. Read with
`SymbolicUtils.getmetadata(sys, OpticalComponent, nothing)`.
"""
struct OpticalComponent end

# -----------------------------------------------------------------------------
# Propagation representations
# -----------------------------------------------------------------------------

"""
    PropagationRep

How `FreeSpace` propagation is represented in the time domain. Passed to
`FreeSpace(rep = ...)` and threaded down from the top of a model, so switching
a whole model between representations is a one-argument change per analysis:

  - [`Algebraic`](@ref) — static phase rotation, the adiabatic limit.
  - [`PadeDelay`](@ref) — order-n all-pass rational approximation of the delay.

For frequency-domain questions that need the *exact* transcendental delay, do
not rationalize at all: lower the model with [`to_scattering`](@ref).
"""
abstract type PropagationRep end

"""
    Algebraic()

Propagation as the static phase rotation `R(2πL/λ)` — no dynamics, no states.
Exact at DC; the right default for audio-band studies of short cavities where
every analysis frequency is far below `c/L`.
"""
struct Algebraic <: PropagationRep end

"""
    PadeDelay(order)

Order-n all-pass rational approximation of the free-space delay `e^{-sτ}`,
`τ = L/c`: a cascade of `order` first-order all-pass sections

    H₁(s) = (1 - sθ/2)/(1 + sθ/2),    θ = τ/order

applied identically to both quadratures, with the static rotation `R(2πL/λ)`
as a separate algebraic factor. `|H₁(iω)| ≡ 1` exactly, so passivity holds at
every order; only phase is approximate (error `~(ωθ)³/12` per section, total
`∝ order⁻²` at fixed `ωτ`). Cavity FSR physics becomes representable once
`order ≳ ωτ` across the band of interest. Costs `2·order` states per
`FreeSpace` instance.
"""
struct PadeDelay <: PropagationRep
    order::Int
    function PadeDelay(order::Int)
        order ≥ 1 || throw(ArgumentError("PadeDelay order must be ≥ 1, got $order"))
        return new(order)
    end
end

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
    meta = [OpticalComponent => (kind = :Mirror, r = sqrt(R), t = sqrt(Tval))]
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

    return System(equations, t, [], pars; name, systems, metadata = meta)
end

# -----------------------------------------------------------------------------
# FreeSpace
# -----------------------------------------------------------------------------
#
# Free-space propagation over length L. Unidirectional — to wire two mirrors
# face-to-face you instantiate one `FreeSpace` for the right-going field and a
# second one for the left-going field.
#
# The physics is `out(t) = R(2πL/λ) · in(t - L/c)`: a static phase rotation of
# the envelope times a pure delay. The rotation is always kept as an exact
# algebraic factor; `rep` chooses how the delay is realized (see
# `PropagationRep`). Because every all-pass section has unit magnitude at all
# frequencies, energy conservation holds for every `rep` at every order.

@component function FreeSpace(; name, L = 1.0, λ = 1064e-9, c = 299792458.0,
        rep::PropagationRep = Algebraic())
    meta = [OpticalComponent =>
        (kind = :FreeSpace, L = float(L), λ = float(λ), c = float(c))]

    pars = @parameters begin
        L = L, [description = "Propagation length [m]"]
        λ = λ, [description = "Carrier wavelength [m]"]
        c = c, [description = "Speed of light [m/s]"]
    end

    systems = @named begin
        port_in = OpticalIn()
        port_out = OpticalOut()
    end

    # Static envelope rotation by φ = k L = (2π/λ) L, applied ahead of the
    # (rep-dependent) delay realization. u holds the rotated input field.
    φ = 2π * L / λ
    u = (cos(φ) * port_in.a[1] - sin(φ) * port_in.a[2],
         sin(φ) * port_in.a[1] + cos(φ) * port_in.a[2])

    if rep isa Algebraic
        equations = Equation[
            port_out.a[1] ~ u[1],
            port_out.a[2] ~ u[2],
        ]
        return System(equations, t, [], pars; name, systems, metadata = meta)
    end

    rep::PadeDelay
    n = rep.order
    θ = L / (c * n)             # per-section delay
    @variables x(t)[1:n, 1:2] [description = "All-pass section states"]

    # Cascade of n first-order all-pass sections per quadrature:
    #     ẋⱼ = -(2/θ)(xⱼ - uⱼ₋₁),   yⱼ = 2xⱼ - uⱼ₋₁
    equations = Equation[]
    for j in 1:n
        push!(equations, D(x[j, 1]) ~ -(2 / θ) * (x[j, 1] - u[1]))
        push!(equations, D(x[j, 2]) ~ -(2 / θ) * (x[j, 2] - u[2]))
        u = (2x[j, 1] - u[1], 2x[j, 2] - u[2])
    end
    push!(equations, port_out.a[1] ~ u[1])
    push!(equations, port_out.a[2] ~ u[2])

    return System(equations, t, [x], pars; name, systems, metadata = meta,
        guesses = [x => zeros(n, 2)])
end

# -----------------------------------------------------------------------------
# Laser (constant-amplitude source)
# -----------------------------------------------------------------------------
#
# Emits a fixed complex envelope `(A_re, A_im)` from its single output port and
# absorbs whatever comes back on its single input port. The laser does not
# respond to back-action.

@component function Laser(; name, A = 1.0, ϕ = 0.0)
    meta = [OpticalComponent => (kind = :Laser,)]
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

    return System(equations, t, [], pars; name, systems, metadata = meta)
end

# -----------------------------------------------------------------------------
# LaserModulated (constant carrier + small-signal quadrature modulation inputs)
# -----------------------------------------------------------------------------
#
# Same as `Laser`, but exposes small-signal perturbations of both carrier-frame
# quadratures as variables marked `[input = true]`: `δA` in phase with the
# carrier (amplitude modulation) and `δQ` in quadrature (phase modulation, to
# first order). The emitted envelope is `(A + δA + i δQ) e^{iϕ}`. Defaulting to
# zero, they do nothing at the DC operating point — but passed as inputs to
# `linearize`, they span the full quadrature pair, so port-to-port transfer
# functions come out as complete 2×2 blocks.

@component function LaserModulated(; name, A = 1.0, ϕ = 0.0)
    meta = [OpticalComponent => (kind = :Laser,)]
    pars = @parameters begin
        A = A, [description = "Carrier field amplitude"]
        ϕ = ϕ, [description = "Carrier field phase [rad]"]
    end

    @variables δA(t) = 0.0 [input = true,
        description = "Small-signal amplitude-quadrature perturbation"]
    @variables δQ(t) = 0.0 [input = true,
        description = "Small-signal phase-quadrature perturbation"]

    systems = @named begin
        port_out = OpticalOut()
        port_in = OpticalIn()  # absorbs reflected light
    end

    equations = Equation[
        port_out.a[1] ~ (A + δA) * cos(ϕ) - δQ * sin(ϕ),
        port_out.a[2] ~ (A + δA) * sin(ϕ) + δQ * cos(ϕ),
    ]

    return System(equations, t, [δA, δQ], pars; name, systems, metadata = meta)
end

# -----------------------------------------------------------------------------
# Vacuum (zero-amplitude source / beam dump)
# -----------------------------------------------------------------------------
#
# Emits zero field on its output port. Used to terminate an unused side of a
# mirror so the network has no dangling directed connectors.

@component function Vacuum(; name)
    meta = [OpticalComponent => (kind = :Vacuum,)]
    systems = @named begin
        port_out = OpticalOut()
        port_in = OpticalIn()
    end

    equations = Equation[
        port_out.a[1] ~ 0,
        port_out.a[2] ~ 0,
    ]

    return System(equations, t, [], []; name, systems, metadata = meta)
end
