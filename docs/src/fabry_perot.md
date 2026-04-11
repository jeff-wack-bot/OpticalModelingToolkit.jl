# The Fabry-Perot cavity

This page documents the optical components defined in `src/components.jl`,
the connectors in `src/connectors.jl`, and the two Fabry-Perot tests in
`test/test_phase1a_fabry_perot.jl` and `test/test_phase1b_fabry_perot_ac.jl`.

## Connectors

```julia
@connector function OpticalIn(;  name, guess = zeros(2)) ... end
@connector function OpticalOut(; name, guess = zeros(2)) ... end
```

Each connector carries a length-2 real array — the real and imaginary parts
of the slowly-varying complex envelope amplitude at the carrier frequency.
They are directed (`[input = true]` and `[output = true]`), modeled on
`RealInputArray`/`RealOutputArray` from `ModelingToolkitStandardLibrary.Blocks`.

A bidirectional optical port is built from one of each. For example, the
front face of a `Mirror` exposes both `f_in::OpticalIn` and `f_out::OpticalOut`.
Two components are wired face-to-face with two `connect` calls — one for
each direction:

```julia
connect(m1.f_out, L_fwd.port_in)
connect(L_fwd.port_out, m2.b_in)

connect(m2.b_out, L_bwd.port_in)
connect(L_bwd.port_out, m1.f_in)
```

**Why directed instead of acausal?** An acausal connector with a flow
variable forces MTK to invent a "directional crossover" inside `connect()`,
which is custom semantics on top of the standard library. Directed
sub-ports get the same Finesse-style port-node directional model "for free"
out of plain `output ~ input` equality, with zero new connector
infrastructure to test.

**Why a length-2 real array instead of a complex scalar?** MTK Standard
Library is real-valued throughout. Splitting `(real, imag)` into two reals
means every component composes with `Mechanical`, `Electrical`, `Blocks`,
etc., without a complex-arithmetic shim.

## Components

### `Mirror(R)`

Two-port lossless dielectric scatterer. Each face has `*_in` and `*_out`
sub-ports. Scattering equations:

```
b_out.a = mt * f_in.a + mr * b_in.a
f_out.a = -mr * f_in.a + mt * b_in.a
```

with `mr = sqrt(R)`, `mt = sqrt(1 − R)`. The scattering matrix
`S = [[mt, mr]; [-mr, mt]]` satisfies `S†S = (mr² + mt²) I`, so the mirror
is unitary when `R + T = 1`. The minus sign is the standard convention
that makes a real-valued lossless dielectric mirror unitary.

The parameter symbols are named `mr`/`mt` rather than `r`/`t` because `t`
would shadow MTK's time variable inside the component constructor.

### `Length(L, λ)`

Algebraic phase rotation: `out.a = exp(ikL) * in.a` with `k = 2π/λ`,
expanded to its real and imaginary components. Unidirectional — to wire
two mirrors face-to-face you instantiate one `Length` for each direction.

### `LengthDelay(L, λ, c)`

Single-pole approximation of free-space delay:

```
τ ȯ + o = exp(ikL) i,    τ = L / c
```

At DC this matches `Length` exactly. At higher modulation frequencies the
response rolls off as `1/(1 + iωτ)`, a one-pole approximation of the true
`e^{−iωτ}` delay.

**Why have both `Length` and `LengthDelay`?** A purely algebraic network
has zero state variables, so `linearize_symbolic` returns only `D` and
`H(iω)` is constant. Promoting propagation to a state-bearing form is what
makes the AC analysis non-trivial; the algebraic version is kept for the
DC test because it has no temporal artifacts to confound the steady-state
comparison. The single-pole form is a deliberate stand-in — good enough to
exercise the linearization machinery and recover an analytical transfer
function, but not a faithful distributed delay model.

### `Laser(A, ϕ)` and `LaserModulated(A, ϕ)`

`Laser` emits a fixed complex envelope `(A cos ϕ, A sin ϕ)` from `port_out`
and absorbs whatever returns on `port_in`. `LaserModulated` is the same
plus a small-signal in-phase amplitude perturbation `δA(t)` declared with
`[input = true]`, defaulting to zero. The variable does nothing at the DC
operating point; its purpose is to be a designated linearization input.

### `Vacuum()`

Emits zero on `port_out`, absorbs whatever arrives on `port_in`. Used to
terminate an unused mirror face so the network has no dangling directed
ports.

## The Fabry-Perot tests

Both tests build the same topology:

```
Laser ──> M1 ──L_fwd──> M2 ──> Vacuum
            ^             |
            └────L_bwd────┘
```

— two mirrors wired face-to-face through two unidirectional length elements
(one per direction), with a `Vacuum` terminator on the far side of M2.

### `test_phase1a_fabry_perot.jl`

Algebraic `Length`, constant `Laser`. The full network is purely
algebraic, so `mtkcompile` produces an equation-only system that is solved
as a `NonlinearProblem` with `NewtonRaphson`. The test computes:

- `m1.f_out.a` — the right-going intracavity field at M1's front face,
  compared against `t1 A / (1 + r1 r2 exp(2ikL))` (the `+` follows from
  our scattering convention)
- `m1.b_out.a` — the field reflected back toward the laser
- `m2.f_out.a` — the field transmitted out the far side of M2
- `|E_refl|² + |E_trans|² ≈ |A|²` — energy conservation

All field comparisons hold to `rtol = 1e-10`.

A small wart: after `mtkcompile`, two scalar array components remain as
free unknowns even though `OpticalIn`/`OpticalOut` declare `guesses = [a
=> zeros(2)]`. The connector-level guess does not propagate to scalar
components after structural simplification, so the test passes
`[u => 0.0 for u in unknowns(sys)]` to `NonlinearProblem` explicitly. Zero
is a valid initial guess for a linear algebraic system.

### `test_phase1b_fabry_perot_ac.jl`

Dynamic `LengthDelay`, `LaserModulated`. The cavity now has 4 state
variables (two quadratures × two length elements). The test calls

```julia
(; A, B, C, D), ssys = ModelingToolkit.linearize_symbolic(
    fp,
    [laser.δA],                           # input
    [m1.f_out.a[1], m1.f_out.a[2]],       # outputs (real, imag quadratures)
)
```

and asserts shape `(4×4, 4×1, 2×4, 2×1)`. The symbolic matrices are
substituted to numbers by walking `parameters(ssys)` and binding each
parameter to its constructor-time `getdefault` value via
`ModelingToolkit.fixpoint_sub`. The system is then verified to be
asymptotically stable: every eigenvalue of `A` has strictly negative real
part.

The numerical transfer function is

```
H_num(ω) = v[1] + i v[2],   v = C (iωI − A)⁻¹ B + D
```

— the two real outputs (real and imaginary quadratures of the field
response) are recombined into a single complex field transfer function.
This is compared against the analytical single-pole-delay cavity response

```
H(ω) = t1 (1 + iωτ)² / [(1 + iωτ)² + r1 r2 exp(2ikL)]
```

at 25 logarithmically-spaced frequencies from 100 Hz to 100 MHz, at
`rtol = 1e-8`. In practice the maximum relative error across the sweep
is ~`5 × 10⁻¹⁴`. The DC value is checked separately against
`t1 / (1 − r1 r2)`.

The cavity length is chosen as an odd multiple of `λ/4` so that
`exp(2ikL) = −1` (resonance under our scattering convention) and the
intracavity buildup is real to within the floating-point representation
of `2π L / λ`.
