# Optical ports.
#
# The carried quantity is the slowly-varying complex envelope amplitude of the
# optical field at the carrier frequency. MTK is real-valued, so a complex
# amplitude travels as a length-2 array `a = [real, imag]` on a *directed*
# connector: `OpticalOut` carries the field leaving a component, `OpticalIn`
# the field entering it. Directed connectors make a connection plain equality
# with a direction, which is exactly a field travelling from one port to the
# next; no new `connect` semantics are needed.
#
# Incoming fields are `irreducible`: they stay unknowns of the compiled system.
# Every component equation then has exactly one variable tearing may solve for,
# its outgoing field, with coefficient 1, and every connection solves `o = i`.
# Without this, tearing solves component equations backwards, dividing by
# coefficients such as `√T` or the `cos 2ϕ`/`sin 2ϕ` of a tuning that can be
# numerically zero (a mirror with `T = 0`, a cavity on resonance) and silently
# make the compiled system singular. The field network
# is solved forward, as in Finesse, never inverted.
#
# A physical port is the pair: `OpticalPort` bundles an `i::OpticalIn` and an
# `o::OpticalOut`, and `link(p, q)` wires two ports face to face with one
# `connect` per direction.

@connector function OpticalIn(; name, guess = zeros(2))
    @variables a(t)[1:2] [
        input = true,
        irreducible = true,
        description = "Complex envelope amplitude (real, imag) entering the component",
    ]
    return System(Equation[], t, [a], []; name, guesses = [a => guess])
end

@connector function OpticalOut(; name, guess = zeros(2))
    @variables a(t)[1:2] [
        output = true,
        description = "Complex envelope amplitude (real, imag) leaving the component",
    ]
    return System(Equation[], t, [a], []; name, guesses = [a => guess])
end

"""
    OpticalPort(; name)

A bidirectional optical port: `i::OpticalIn` carries the field arriving at the
component, `o::OpticalOut` the field leaving it. Every component names its ports
`p1`, `p2`, …; wire two of them together with [`link`](@ref).
"""
@component function OpticalPort(; name)
    systems = @named begin
        i = OpticalIn()
        o = OpticalOut()
    end
    return System(Equation[], t, [], []; name, systems)
end

"""
    link(p, q) -> Vector{Equation}

Wire optical ports `p` and `q` face to face: the field leaving `p` enters `q`,
and the field leaving `q` enters `p`.

```julia
eqs = [link(laser.p1, m1.p1); link(m1.p2, space.p1); link(space.p2, m2.p1)]
```
"""
link(p, q) = [connect(p.o, q.i), connect(q.o, p.i)]

"""
    expose(outer, inner) -> Vector{Equation}

Inside a subassembly, make the port `inner` of one of its components available
as the subassembly's own port `outer`: the field arriving at `outer` passes in
to `inner`, and the field leaving `inner` passes out through `outer`.

```julia
@component function Cavity(; name, L)
    @named m1 = Mirror(); @named space = FreeSpace(; L); @named m2 = Mirror()
    @named p1 = OpticalPort(); @named p2 = OpticalPort()
    eqs = [expose(p1, m1.p1); link(m1.p2, space.p1); link(space.p2, m2.p1); expose(p2, m2.p2)]
    return System(eqs, t; name, systems = [m1, space, m2, p1, p2])
end
```
"""
expose(outer, inner) = [connect(outer.i, inner.i), connect(inner.o, outer.o)]

# Components are written in complex arithmetic and realified here: `field`
# reads a connector as a `Complex{Num}`, `emit` sets a connector from one.

field(c) = Complex(c.a[1], c.a[2])
emit(c, z) = [c.a[1] ~ real(z), c.a[2] ~ imag(z)]
