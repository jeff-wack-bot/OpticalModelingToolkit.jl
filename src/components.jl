# Optical components.
#
# Plane-wave (TEM₀₀) optics on the complex envelope. Every component is written
# once, in complex arithmetic, as a scattering table: a list of
# `(out_port, in_port) => amplitude` entries saying how much of the field
# entering `in_port` leaves through `out_port`. The same table function is
# evaluated twice:
#
#   - with symbolic parameters, giving the component's MTK equations (realified
#     by `emit`, so every complex product becomes a 2×2 real block), and
#   - with numbers, giving the `OpticalComponent` metadata that `to_scattering`
#     reads to lower a composed model to the exact frequency-domain form.
#
# Free-space propagation is the one place dynamics enter, and *how* the delay is
# represented is an analysis-time choice threaded in through
# `FreeSpace(rep = ...)`, never part of the component's identity.
#
# Sign convention: a mirror reflects with `+r` on its `p1` face and `-r` on its
# `p2` face, and transmits with `+t` both ways. The scattering matrix is real
# and orthogonal for `R + T = 1`. Propagation multiplies by `e^{+ikL}`.

# -----------------------------------------------------------------------------
# Component metadata
# -----------------------------------------------------------------------------

"""
    OpticalComponent

`System`-metadata key under which every optical component records its role in
the scattering network. `to_scattering` reads it to lower a composed model; the
value is one of

  - [`Scatterer`](@ref): a static linear map between ports,
  - [`Propagator`](@ref): a two-port delay with a phase rotation,
  - [`Terminal`](@ref): a one-port source or sink, which becomes an external
    channel.

A new component joins the frequency-domain lowering by setting this metadata.
"""
struct OpticalComponent end

"""
    Scatterer(S)

Role of a component whose ports are related by the static complex scattering
table `S`, a vector of `(out_port, in_port, amplitude)` triples.
"""
struct Scatterer
    S::Vector{Tuple{Symbol, Symbol, ComplexF64}}
end
Scatterer(S::AbstractVector{<:Pair}) = Scatterer([(o, i, ComplexF64(c)) for ((o, i), c) in S])

"""
    Propagator(τ, φ)

Role of a two-port (`p1`, `p2`) that delays the field by `τ` and rotates it by
`e^{iφ}` in both directions.
"""
struct Propagator
    τ::Float64
    φ::Float64
end

"""
    Terminal()

Role of a one-port (`p1`) source or sink. In the scattering form, the field it
emits is an external input channel and the field it receives is an external
output channel, both named after the component.
"""
struct Terminal end

# `n` optical ports named p1, …, pn, as a NamedTuple.
function optical_ports(n)
    names = ntuple(k -> Symbol(:p, k), n)
    return NamedTuple{names}(map(name -> OpticalPort(; name), names))
end

# Realified equations `p.o = Σ S[p, q] · q.i` for every port `p`, from a table
# of `(out, in) => amplitude` pairs over the `ports` NamedTuple.
function scattering_equations(ports::NamedTuple, S)
    eqs = Equation[]
    for (name, port) in pairs(ports)
        z = zero(Complex{Num})
        for ((o, i), c) in S
            o == name && (z += c * field(ports[i].i))
        end
        append!(eqs, emit(port.o, z))
    end
    return eqs
end

# -----------------------------------------------------------------------------
# Propagation representations
# -----------------------------------------------------------------------------

"""
    PropagationRep

How `FreeSpace` propagation is represented in the time domain. Passed to
`FreeSpace(rep = ...)` and threaded down from the top of a model, so switching
a whole model between representations is a one-argument change per analysis:

  - [`Algebraic`](@ref): static phase rotation, the adiabatic limit.
  - [`PadeDelay`](@ref): order-n all-pass rational approximation of the delay.

For frequency-domain questions that need the *exact* transcendental delay, do
not rationalize at all: lower the model with [`to_scattering`](@ref).
"""
abstract type PropagationRep end

"""
    Algebraic()

Propagation as the static phase rotation `e^{ikL}`, with no dynamics and no
states. Exact at DC; the right default for audio-band studies of short
cavities where every analysis frequency is far below `c/L`.
"""
struct Algebraic <: PropagationRep end

"""
    PadeDelay(order)

Order-n all-pass rational approximation of the free-space delay `e^{-sτ}`,
`τ = L/c`: a cascade of `order` first-order all-pass sections

    H₁(s) = (1 - sθ/2)/(1 + sθ/2),    θ = τ/order

applied to the complex envelope after the static rotation `e^{ikL}`.
`|H₁(iω)| ≡ 1` exactly, so passivity holds at every order; only the phase is
approximate (error `~(ωθ)³/12` per section, total `∝ order⁻²` at fixed `ωτ`).
Cavity FSR physics becomes representable once `order ≳ ωτ` across the band of
interest. Costs `4·order` real states per `FreeSpace` (two quadratures, two
directions).
"""
struct PadeDelay <: PropagationRep
    order::Int
    function PadeDelay(order::Int)
        order ≥ 1 || throw(ArgumentError("PadeDelay order must be ≥ 1, got $order"))
        return new(order)
    end
end

# Realify an all-pass cascade on complex input `u`, with states `x[1:n, 1:2]`:
#     ẋⱼ = -(2/θ)(xⱼ - uⱼ₋₁),   uⱼ = 2xⱼ - uⱼ₋₁
function allpass_cascade(u, x, n, θ)
    eqs = Equation[]
    for j in 1:n
        xj = Complex(x[j, 1], x[j, 2])
        ẋ = -(2 / θ) * (xj - u)
        push!(eqs, D(x[j, 1]) ~ real(ẋ), D(x[j, 2]) ~ imag(ẋ))
        u = 2xj - u
    end
    return eqs, u
end

# -----------------------------------------------------------------------------
# Scatterers
# -----------------------------------------------------------------------------

mirror_scattering(r, τ, θ = 0) = [
    (:p1, :p1) => r * cis(θ), (:p1, :p2) => τ,
    (:p2, :p2) => -r * cis(-θ), (:p2, :p1) => τ,
]

"""
    Mirror(; R = 0.99, T = 1 - R)

Two-port mirror with power reflectivity `R` and transmissivity `T`, lossy when
`R + T < 1`. Reflection is `+√R` on the `p1` face and `-√R` on the `p2` face;
transmission is `√T` both ways.
"""
@component function Mirror(; name, R = 0.99, T = 1 - R)
    meta = [OpticalComponent => Scatterer(mirror_scattering(sqrt(R), sqrt(T)))]
    pars = @parameters begin
        R = R, [description = "Power reflectivity"]
        T = T, [description = "Power transmissivity"]
    end
    ports = optical_ports(2)
    eqs = scattering_equations(ports, mirror_scattering(sqrt(R), sqrt(T)))
    return System(eqs, t, [], pars; name, systems = collect(ports), metadata = meta)
end

"""
    MovableMirror(; R = 0.99, T = 1 - R, λ = 1064e-9, c = 299792458.0)

A [`Mirror`](@ref) with a microscopic longitudinal position `x` (an
`[input = true]` variable, positive toward the `p2` side) and the radiation
pressure force `F` it feels (an `[output = true]` variable, positive toward the
`p2` side). Bind `x` with an equation, e.g. to a mechanical model driven by
`F`, or name it as a linearization input.

Moving by `x` lengthens the `p1`-side round trip and shortens the `p2`-side one:
the `p1` reflection picks up `e^{+2ikx}` and the `p2` reflection `e^{-2ikx}`.
The full nonlinear dependence is kept, so a length-control loop can sweep
through fringes. With field amplitudes in `√W`, the force is

    F = (|p1.i|² + |p1.o|² - |p2.i|² - |p2.o|²) / c.

At `x = 0` the scattering is exactly that of `Mirror`.
"""
@component function MovableMirror(; name, R = 0.99, T = 1 - R, λ = 1064e-9,
        c = 299792458.0)
    meta = [OpticalComponent => Scatterer(mirror_scattering(sqrt(R), sqrt(T)))]
    pars = @parameters begin
        R = R, [description = "Power reflectivity"]
        T = T, [description = "Power transmissivity"]
        λ = λ, [description = "Carrier wavelength [m]"]
        c = c, [description = "Speed of light [m/s]"]
    end
    # guesses, not defaults: x is set by whatever the caller binds it to and F
    # is algebraically determined, so defaults would overdetermine
    # initialization
    vars = @variables begin
        x(t), [input = true, description = "Position [m], positive toward p2"]
        F(t), [output = true, description = "Radiation pressure force [N], positive toward p2"]
    end
    ports = optical_ports(2)
    (; p1, p2) = ports
    eqs = scattering_equations(ports, mirror_scattering(sqrt(R), sqrt(T), 2 * (2π / λ) * x))
    push!(eqs, F ~ (abs2(field(p1.i)) + abs2(field(p1.o)) -
                    abs2(field(p2.i)) - abs2(field(p2.o))) / c)
    return System(eqs, t, vars, pars; name, systems = collect(ports), metadata = meta,
        guesses = [x => 0.0, F => 0.0])
end

beamsplitter_scattering(r, τ) = [
    (:p1, :p2) => r, (:p2, :p1) => r, (:p3, :p4) => -r, (:p4, :p3) => -r,
    (:p1, :p3) => τ, (:p3, :p1) => τ, (:p2, :p4) => τ, (:p4, :p2) => τ,
]

"""
    Beamsplitter(; R = 0.5, T = 1 - R)

Four-port beamsplitter. `p1` and `p2` are on the front face and reflect into
each other with `+√R`; `p3` and `p4` are on the back face and reflect into each
other with `-√R`. Transmission (`√T`) couples `p1 ↔ p3` and `p2 ↔ p4`.
"""
@component function Beamsplitter(; name, R = 0.5, T = 1 - R)
    meta = [OpticalComponent => Scatterer(beamsplitter_scattering(sqrt(R), sqrt(T)))]
    pars = @parameters begin
        R = R, [description = "Power reflectivity"]
        T = T, [description = "Power transmissivity"]
    end
    ports = optical_ports(4)
    eqs = scattering_equations(ports, beamsplitter_scattering(sqrt(R), sqrt(T)))
    return System(eqs, t, [], pars; name, systems = collect(ports), metadata = meta)
end

# -----------------------------------------------------------------------------
# FreeSpace
# -----------------------------------------------------------------------------

"""
    FreeSpace(; L = 1.0, λ = 1064e-9, c = 299792458.0, rep = Algebraic())

Free-space propagation over length `L` between ports `p1` and `p2`, in both
directions: `out(t) = e^{ikL} in(t - L/c)` with `k = 2π/λ`. The rotation is kept
exact; `rep` (a [`PropagationRep`](@ref)) chooses how the delay is realized.
"""
@component function FreeSpace(; name, L = 1.0, λ = 1064e-9, c = 299792458.0,
        rep::PropagationRep = Algebraic())
    meta = [OpticalComponent => Propagator(L / c, 2π * L / λ)]
    pars = @parameters begin
        L = L, [description = "Propagation length [m]"]
        λ = λ, [description = "Carrier wavelength [m]"]
        c = c, [description = "Speed of light [m/s]"]
    end
    (; p1, p2) = optical_ports(2)
    P = cis(2π * L / λ)
    u12, u21 = P * field(p1.i), P * field(p2.i)

    if rep isa Algebraic
        eqs = [emit(p2.o, u12); emit(p1.o, u21)]
        return System(eqs, t, [], pars; name, systems = [p1, p2], metadata = meta)
    end

    n = rep.order
    θ = L / (c * n)
    @variables begin
        x12(t)[1:n, 1:2], [description = "All-pass states, p1 → p2"]
        x21(t)[1:n, 1:2], [description = "All-pass states, p2 → p1"]
    end
    eqs12, y12 = allpass_cascade(u12, x12, n, θ)
    eqs21, y21 = allpass_cascade(u21, x21, n, θ)
    eqs = [eqs12; eqs21; emit(p2.o, y12); emit(p1.o, y21)]
    return System(eqs, t, [x12, x21], pars; name, systems = [p1, p2], metadata = meta,
        guesses = [x12 => zeros(n, 2), x21 => zeros(n, 2)])
end

# -----------------------------------------------------------------------------
# Terminals: one-port sources and sinks
# -----------------------------------------------------------------------------

"""
    Laser(; A = 1.0, ϕ = 0.0)

Emits the constant envelope `A e^{iϕ}` from `p1` and absorbs whatever arrives.
"""
@component function Laser(; name, A = 1.0, ϕ = 0.0)
    pars = @parameters begin
        A = A, [description = "Field amplitude [√W]"]
        ϕ = ϕ, [description = "Field phase [rad]"]
    end
    p1 = OpticalPort(; name = :p1)
    return System(emit(p1.o, A * cis(ϕ)), t, [], pars; name, systems = [p1],
        metadata = [OpticalComponent => Terminal()])
end

"""
    LaserModulated(; A = 1.0, ϕ = 0.0)

A [`Laser`](@ref) whose envelope `(A + δA + i δQ) e^{iϕ}` carries small-signal
modulation inputs: `δA` in phase with the carrier (amplitude modulation) and
`δQ` in quadrature (phase modulation, to first order). Both are
`[input = true]` variables: bind them with equations or name them as
linearization inputs, where together they span the full quadrature pair.
"""
@component function LaserModulated(; name, A = 1.0, ϕ = 0.0)
    pars = @parameters begin
        A = A, [description = "Carrier amplitude [√W]"]
        ϕ = ϕ, [description = "Carrier phase [rad]"]
    end
    vars = @variables begin
        δA(t) = 0.0, [input = true, description = "Amplitude-quadrature modulation"]
        δQ(t) = 0.0, [input = true, description = "Phase-quadrature modulation"]
    end
    p1 = OpticalPort(; name = :p1)
    return System(emit(p1.o, (A + δA + im * δQ) * cis(ϕ)), t, vars, pars; name,
        systems = [p1], metadata = [OpticalComponent => Terminal()])
end

"""
    Vacuum()

Emits zero field from `p1` and absorbs whatever arrives: terminates an unused
port so the network has no dangling connectors.
"""
@component function Vacuum(; name)
    p1 = OpticalPort(; name = :p1)
    return System(emit(p1.o, 0), t, [], []; name, systems = [p1],
        metadata = [OpticalComponent => Terminal()])
end

"""
    Photodetector()

Perfect absorber on `p1` that reports the arriving power `P = |p1.i|²` as an
`[output = true]` variable.
"""
@component function Photodetector(; name)
    # a guess, not a default: P is algebraically determined by the field
    @variables P(t) [output = true, description = "Detected power [W]"]
    p1 = OpticalPort(; name = :p1)
    eqs = [P ~ abs2(field(p1.i)); emit(p1.o, 0)]
    return System(eqs, t, [P], []; name, systems = [p1],
        metadata = [OpticalComponent => Terminal()], guesses = [P => 0.0])
end
