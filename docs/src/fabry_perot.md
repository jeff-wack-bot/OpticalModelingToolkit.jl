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

### `FreeSpace(L, λ, c, rep)`

Free-space propagation over length `L`: physically
`out(t) = R(2πL/λ) · in(t − L/c)` — a static phase rotation of the envelope
times a pure delay. Unidirectional — to wire two mirrors face-to-face you
instantiate one `FreeSpace` for each direction.

The rotation is always exact and algebraic; `rep::PropagationRep` chooses how
the delay is represented, per analysis, threaded in from the top of a model
(how propagation is represented is an analysis-time choice, not an element
identity — `plans/PLAN_rationalize_late.md` §2):

- `Algebraic()` (default) — no dynamics, no states; the adiabatic limit,
  exact at DC.
- `PadeDelay(n)` — a cascade of `n` first-order **all-pass** sections
  `H₁(s) = (1 − sθ/2)/(1 + sθ/2)`, `θ = τ/n`, applied identically to both
  quadratures. `|H₁(iω)| ≡ 1` exactly, so passivity/energy conservation holds
  at every order; only phase is approximate (per-section error `(ωθ)³/12`,
  total `∝ n⁻²` at fixed `ωτ`). Cavity FSR physics is representable once
  `n ≳ ωτ` over the band of interest — a convergence knob, not a wall. Costs
  `2n` states per instance.

For exact delays in the frequency domain, skip rationalization entirely:
[`to_scattering`](@ref) lowers the composed model to the transcendental
scattering form (each `FreeSpace` an exact `e^{-sτ}`), and with Etalon loaded
`Etalon.ScatteringSystem(model)` evaluates it directly.

### `Laser(A, ϕ)` and `LaserModulated(A, ϕ)`

`Laser` emits a fixed complex envelope `(A cos ϕ, A sin ϕ)` from `port_out`
and absorbs whatever returns on `port_in`. `LaserModulated` is the same plus
small-signal perturbations of both carrier-frame quadratures declared with
`[input = true]` and defaulting to zero: `δA(t)` in phase with the carrier
(amplitude modulation) and `δQ(t)` in quadrature (phase modulation to first
order), emitting `(A + δA + i δQ) e^{iϕ}`. They do nothing at the DC
operating point; their purpose is to be designated linearization inputs that
span the full quadrature pair.

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

`FreeSpace` in its algebraic representation, constant `Laser`. The full
network is purely algebraic, so `mtkcompile` produces an equation-only system that is solved
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

`FreeSpace(rep = PadeDelay(1))`, `LaserModulated`. The cavity has 4 state
variables (two quadratures × two propagation directions × one all-pass
section). The AC path is the ecosystem route documented in
[Solver strategy](solvers.md):

```julia
mats, _ = linearize_proper(fp, [laser.δA, laser.δQ],
    [m1.f_out.a[1], m1.f_out.a[2]])
sys = ss(mats...)
freqresp(sys, ωs)
```

The two real outputs (quadratures of the field response to the in-phase
drive) are recombined into a single complex field transfer function and
compared against the representation-parameterized analytic reference

```
H(s) = t1 / (1 + r1 r2 P(s)²),    P(s) = e^{ikL} (1 − sτ/2)/(1 + sτ/2)
```

at 25 logarithmically-spaced frequencies from 100 Hz to 100 MHz. Because the
all-pass section has exact unit magnitude, `|H|` is exact at every frequency
and only the phase carries the representation error. `poles` and `dcgain`
(`t1/(1 − r1 r2)` on resonance) are checked through ControlSystemsBase, and
the whole sweep is cross-checked against the hand-rolled symbolic oracle
(`linearize_symbolic` → `fixpoint_sub` → `C(iωI−A)⁻¹B + D`), which is kept
for pedagogy only.

The cavity length is chosen as an odd multiple of `λ/4` so that
`exp(2ikL) = −1` (resonance under our scattering convention) and the
intracavity buildup is real to within the floating-point representation
of `2π L / λ`.

Beyond these two, `test/test_propagation.jl` carries the representation
validation ladder (all-pass unit magnitude, `n⁻²` phase-error scaling,
Fabry-Perot convergence to the exact delay response, FSR-peak location,
audio-band representation invariance) and `test/test_scattering.jl` the
transcendental bridge and the A-matrix two-photon keystone test.
