# Shared models with Finesse.jl

A model of an interferometer is written once, in `models/<name>.kat`, and run by
two simulators:

- [Finesse.jl](https://git.mccullerlab.com/software/finesse.jl), a port of
  FINESSE 3, solves its frequency response;
- OMT builds the same network as a ModelingToolkit model with
  [`build_kat`](@ref), and simulates it in the time domain: optics, suspension
  and electronics in one DAE.

The two are complementary. Finesse linearizes about a working operating point,
so it answers questions about small motions around a lock: transfer functions,
sensing matrices, noise. Everything else is a time-domain question: acquiring a
lock, losing it, a mirror swinging through resonance, a saturating actuator. The
examples page *Shared models: one file, two simulators* shows both on the same
files, and checks that they agree.

## Why KatScript is the single source of truth

The model has to be readable by both codes, so there are three choices: a
Julia netlist that writes KatScript for Finesse, KatScript that OMT reads, or a
third format both translate from. KatScript wins because it is already the
language interferometer models are written in: the LIGO and Virgo models,
Finesse's documentation and its users all speak it, and Finesse.jl parses all of
it, so a model is never limited by what OMT supports. OMT reads a subset.

The subset is deliberately strict. A line OMT does not understand is an error,
not a guess: a `link(...)` that OMT silently skipped would drop a servo
connection and produce a model that runs and is wrong. Analysis lines (`fsig`,
`frequency_response`, `xaxis`) are the exception; they say what Finesse should
compute, not what the network is, so OMT skips them.

| KatScript | OMT |
|:----------|:----|
| `laser l P=1` | `Laser(A = √P)` |
| `mirror m R= T= L= phi=` | `Mirror`, or `MovableMirror` if a pendulum hangs it |
| `space s a.pN b.pN L=` | `FreeSpace(L, rep)`, linked to both ports |
| `pendulum p m.mech mass= fz= Qz=` | ``M\ddot z + M\omega_0 \dot z/Q + M\omega_0^2 (z - z_0) = F`` |
| `readout_dc r m.pN.o` | ``P = \lvert E\rvert^2 - P_0`` at that node |
| `zpk f [zeros] [poles] gain=` | [`ZPK`](@ref), Finesse's normalization |
| `wire w a b gain=`, `link(a, b)` | the signal at `a`, times `gain`, added into `b` |

Numbers must be plain literals; expressions, units and references are Finesse
features OMT does not evaluate, and are errors.

## Conventions

The two codes describe the same physics in different coordinates. Each mapping
below was established by a cross-check in `test/finesse/`, not assumed:

- **A Finesse mirror is a turned-around OMT mirror.** Finesse reflects with
  ``-r e^{2i\phi}`` on `p1`, OMT with ``+r e^{2i\phi}``. So Finesse `pN` is OMT
  `p(3-N)`, and Finesse's tuning `phi` [deg] is OMT's ``\phi`` [rad] with the
  opposite sign.
- **A suspended mirror's tuning follows ``-z``.** Found by the detuned suspended
  cavity, where the optical spring is nonzero; at resonance the sign is invisible.
- **Open ports are dark.** Finesse leaves an unconnected node without light; OMT
  terminates it with a `Vacuum`.
- **The operating point is Finesse's.** Finesse linearizes with every mirror at
  its tuning, so the static radiation force never moves a suspended mirror, and a
  readout's DC power is not part of its signal. OMT puts its model in that state:
  [`build_kat`](@ref) solves the static network once, sets the suspension's rest
  position so the radiation force is balanced at ``z = 0``, and subtracts each
  readout's DC power. Error signals and filter outputs are then zero at the
  operating point, so a time-domain simulation started there stays there.

## Keeping the references current

Finesse.jl is not registered, and it lives on a server CI cannot reach. So its
results are stored: `models/reference/<name>.txt` holds Finesse's responses of
`models/<name>.kat` at 31 frequencies, and the main test suite and these docs
compare OMT against them.

After changing a model:

```
julia --project=test/finesse -e 'using Pkg; Pkg.instantiate()'
julia --project=test/finesse test/finesse/reference.jl     # rewrite the references
julia --project=test/finesse test/finesse/runtests.jl      # check OMT against live Finesse
```

The cross-check suite fails if a stored reference differs from what Finesse
computes now, so a stale reference cannot pass unnoticed.

## Choosing a delay representation

`build_kat(text; rep)` sets the representation of every `space`. The choice is
physics, not only speed:

- `Algebraic()` has no delay states: the cavity responds instantly. That is
  exact for the DC operating point and for the *stiffness* of an optical spring.
  It misses the spring's *damping*, which comes from the cavity's storage time.
  On the detuned suspended cavity that is an error of up to ``4\times10^{-5}`` near
  the resonance.
- `PadeDelay(n)` keeps `n` all-pass sections per direction. The cavity's one-way
  delay is 0.33 ns, so even one section's phase error at 100 Hz is
  ``(\omega\tau)^3/12 \sim 10^{-22}``: exact, for these models, at any order.

Delay states make the model very stiff: the cavity's modes sit at
``10^9``–``10^{12}`` /s and the suspension at ``10^2``. Evaluating its frequency
response in Float64, with `ControlSystemsBase.freqresp` or anything else, loses
up to ten digits; on the locked cavity it put a 0.8 % error on the servo filter,
which has nothing to do with the delays. The linearized matrices are exact, so
[`frequency_response`](@ref) evaluates them in `BigFloat` and recovers the full
answer. Use it for any model with `PadeDelay`. A time-domain integration is not
affected; it only has to resolve the stiffness, which a Rosenbrock method does.

## What each code does

| | Finesse.jl | OMT |
|:--|:--|:--|
| Plane-wave carrier fields, DC operating point | yes | yes |
| Higher-order spatial modes, beam tracing, mirror maps | yes | no |
| Linear frequency response, exact delays | yes | yes: `to_scattering` (optics only) and `linearize_proper` (everything) |
| Time domain, nonlinear, transients | no | yes |
| Radiation pressure, optical springs | yes | yes |
| Suspensions | free mass, pendulum, ZPK | any equations: multi-stage, nonlinear |
| Electronics | ZPK, Butterworth, Chebyshev, amplifiers, locks | `ZPK`, and any MTK block: saturation, nonlinearity, logic |
| RF modulation and demodulation | yes | not yet (amplitude and phase inputs only) |
| Quantum noise | yes | no |
| Model language | KatScript, Python-compatible API | Julia, and the KatScript subset above |

## Speed

| step | time |
|:-----|-----:|
| Finesse.jl: frequency response, 31 points | 12.4 ms |
| Finesse.jl: frequency response, 1001 points | 89.5 ms |
| OMT: `build_kat` + `linearize_proper` (compiles the model) | 185.0 ms |
| OMT: `frequency_response` (BigFloat), 31 points | 35.0 ms |
| OMT: `frequency_response` (BigFloat), 1001 points | 1.18 s |
| OMT: `ControlSystemsBase.freqresp` (Float64, inexact here), 1001 points | 1.28 ms |
| OMT: 0.3 s time-domain simulation of the loop | 56.9 ms |

Measured on a Raspberry Pi 5 with Julia 1.13 on `models/locked_cavity.kat`,
by `test/finesse/benchmark.jl`. Each time is a repeat run, so Julia's one-time
compilation of the packages (on the order of a minute here) is excluded.

The two codes do different work. Finesse solves the sparse linear signal problem
directly at each frequency, which is fast and scales gently with the number of
points. OMT builds a symbolic model and linearizes it once, then answers
questions of the result: evaluating its frequency response in extended precision
costs about a millisecond per frequency, and Float64 is a thousand times faster
but wrong for a stiff model like this one. What OMT has that Finesse does not is
the last row: the whole loop, simulated in time, in a fraction of a second.

## Bugs found on the way

Writing the locked cavity found two bugs in Finesse.jl's KatScript front end,
both of which made a servo that drives a mirror's force impossible to write:

- `wire act servo.p2.o m2.mech.F_z`: the wire factory reduced each end to its
  port, so the force node became the six-node `m2.mech` and was refused as
  ambiguous;
- `link(servo.p2.o, m2.mech.F_z)`: `link` passed a space's `L` to the wire it
  built, a `MethodError`.

Each is fixed, with its own regression test, on its own Finesse.jl branch
(`fix/wire-node-ends` and `fix/link-signal-kwargs`), submitted upstream. Until
the `wire` fix is merged, the locked-cavity cross-check needs a Finesse.jl with
it; the stored reference was written with both.
