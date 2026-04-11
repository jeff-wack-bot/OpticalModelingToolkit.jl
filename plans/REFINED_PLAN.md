# OpticalModelingToolkit.jl — Refined Design Plan

This document supersedes the architecture sections of `PROJECT_PLAN.md`. It is the result of studying the source of Finesse 3, ModelingToolkit.jl, and ModelingToolkitStandardLibrary.jl in `context/`, and reflects design decisions made during that review.

## 1. The Architectural Seam

Finesse is a frequency-domain linear field solver with no time-domain integrator. ModelingToolkit is a symbolic DAE compiler with no first-class frequency-domain representation. The two paradigms intersect cleanly at exactly one place:

> **MTK's `linearize_symbolic`** (`ModelingToolkit.jl/src/linearization.jl:544`) extracts a symbolic `(A, B, C, D)` state-space form from a nonlinear DAE around an operating point.

This is the seam OpticalModelingToolkit is designed around. The package uses MTK's DAE machinery to find fixed operating points, then uses `linearize_symbolic` at those points to produce a state-space model from which AC analysis is performed by direct evaluation of `H(iω) = C(iωI − A)⁻¹B + D`. Finesse is used to validate both stages but is not part of the runtime.

## 2. The Hybrid Approach

OpticalModelingToolkit components produce DAE equations in time. The same component definition serves two analysis modes:

- **DC operating point**: solve the DAE for steady state via MTK's nonlinear solver path.
- **Small-signal AC analysis**: call `linearize_symbolic` at the operating point and evaluate `H(iω)` at the frequencies of interest.

## 3. Domain Coverage

The constraint that everything must be expressible as a real-valued MTK DAE determines what OpticalModelingToolkit can and cannot model. The relationship to Finesse is roughly orthogonal — there is large mutual coverage of "control-systems optics," and large parts of each tool that the other cannot reach.

### What OpticalModelingToolkit can model

- Passive linear optical networks (mirrors, beamsplitters, spaces, lenses as thin elements).
- Lasers as carrier-amplitude sources with amplitude/phase modulation inputs.
- Photodetectors as non-loading sensors producing real-valued power signals.
- Mechanical suspensions, free masses, and pendulums via the existing `Mechanical.Translational` and `Mechanical.Rotational` components from MTK Standard Library.
- Radiation-pressure coupling between optical power and mechanical force. 
- Analog electronics (op-amps, filters, drivers) via MTK Standard Library `Electrical` and `Blocks`.
- Discrete-time digital control via MTK's `ShiftIndex(t)` discrete equations, mixed with continuous dynamics in the same model.
- Closed-loop length and frequency control of cavities, including the full electro-mechanical loop.
- Audio-band sideband transfer functions via linearization at the DC operating point.

### What OpticalModelingToolkit cannot (and will not) model

- Higher-order Hermite-Gaussian spatial modes (Phase 5+ research, see §10).
- Gaussian beam tracing, beam parameter propagation, ABCD-matrix mode matching.
- Mode-mismatch loss, alignment sensing, knife-edge scattering, clipping losses — anything that requires a spatial profile.
- Nonlinear optics (second-harmonic generation, parametric down-conversion, Kerr effects).
- Quantum noise as a primitive (radiation-pressure quantum noise, shot noise) — this requires a separate Heisenberg-Langevin layer on top of linearized models, which is a research project.
- Fast carrier dynamics — OpticalModelingToolkit works on the slowly-varying envelope. Anything faster than the simulation timestep is out of reach by construction.

### Where OpticalModelingToolkit exceeds Finesse

These are the things OpticalModelingToolkit can do that Finesse cannot, and which justify building a new tool rather than just driving Finesse from Julia:

- **Time-domain transients with full nonlinear electronics in the loop.** Finesse has no time-domain solver.
- **Digital control loops with arbitrary sample-and-hold and quantization.** Finesse's `Lock` action is a thin abstraction; MTK can model the full DSP chain.
- **Symbolic preprocessing and code generation** for fast repeated AC analysis under parameter sweeps.
- **Hierarchical, parameterizable subassemblies** through MTK's standard composition machinery — a `FabryPerot(L=…, R1=…, R2=…)` module reusable across projects.

### The overlap with Finesse

For passive optical networks driven at DC with audio-band sidebands and mechanical/electrical coupling, the two tools should give bit-identical results within numerical tolerance. **That overlap is the entire validation surface** — every OpticalModelingToolkit component will be tested by building the equivalent Finesse model in Python and comparing.

## 4. Connector Design

Optical ports are built from MTK's directed (signal-flow) connectors. MTK Standard Library `Blocks` provides `RealInput` and `RealOutput` connectors (`MTKStdLib/src/Blocks/utils.jl`) using `[input = true]` and `[output = true]` variable metadata. These are causal: an output connects to an input, and the connection means literal equality with a defined direction. Directionality is enforced by metadata, not by connection assembly.

Each Finesse-style optical port becomes **two directed sub-ports**, an "out" carrying the field leaving the component and an "in" carrying the field entering. Two components are wired face-to-face with two explicit `connect` calls — one for each direction:

```julia
connect(m1.p2_out, m2.p1_in)
connect(m2.p1_out, m1.p2_in)
```

This matches Finesse's directional node model exactly (`finesse3/src/finesse/components/node.py:27-180`, `NodeDirection.INPUT`/`OUTPUT`), uses MTK infrastructure that already exists and is well-tested, and avoids any need for custom directional-crossover semantics in `connect()`.

The Phase 1 spike validates this design end-to-end: build a Beamsplitter with directed ports, run it through `mtkcompile`, and verify the equations look right before generalizing.

## 5. Cross-Domain Components

MTK Standard Library's `EMF` component (`MTKStdLib/src/Electrical/Analog/ideal_components.jl:294-344`) is the canonical example of a multi-domain component: it has both an electrical `OnePort` and a rotational `Flange + Support`, and the constitutive equations couple them (`k * w ~ v`, `flange.tau ~ -k * i`).

Optical components with mechanical or electrical ports follow the same pattern. A `Mirror` with longitudinal radiation-pressure coupling will declare:

- A directed optical port group `(p1_in, p1_out, p2_in, p2_out)`.
- A `Translational.MechanicalPort` `mech` for longitudinal motion.
- Constitutive equations giving the field scattering and the radiation pressure force as a function of incident power.

This is the strongest single argument that OpticalModelingToolkit fits cleanly into the MTK ecosystem rather than fighting it. The cross-domain pattern is already established and tested in the standard library; we are just adding optics as another domain.

## 6. Exporter Strategy

Finesse models are built from Julia by walking the OpticalModelingToolkit component tree and constructing a live Python `finesse.Model` object via PythonCall. Each OpticalModelingToolkit component implements one method:

```julia
to_finesse(component, model::Py)
```

which mutates the Python model in place — adding the component, setting its parameters, declaring its connections.

Human-readable KatScript output is **not** something we need to implement. Finesse already provides `Model.unparse()` (`finesse3/src/finesse/model.py:2552-2566`) which round-trips a `finesse.Model` back to KatScript text. After `to_finesse(...)` builds the model, calling `pyfinesse_model.unparse()` yields a debuggable script for free. This is enough for test fixtures and for humans inspecting what was built.

The Finesse 3 API gotchas already documented in `failing_tests/test_finesse_interop.jl` (Space requires `portA`/`portB`, R+T+L floating-point sensitivity, Laser wavelength kwarg behavior) are the reference for what the exporter must handle correctly.

### Hierarchical name flattening

Finesse names are flat — there is no concept of `cavity.m1`. MTK is hierarchical with `.` separators. The exporter must flatten hierarchical names into Finesse-legal flat names without collisions. The convention will be `_` as the separator (`cavity_m1`), with collision detection during export. This is a small but real piece of work that must live in the exporter from day one.

## 7. Validation Strategy

Two test classes, two tolerance regimes:

- **DC operating point**: `rtol = 1e-8`. Steady-state field amplitudes at every node should agree to near machine precision when both sides use the same physics. Any larger discrepancy is a bug in either the component definition or the exporter.
- **AC transfer function** (Bode magnitude and phase, sweeping a sideband frequency): `rtol = 1e-3`. Both sides have small numerical differences in linearization around the operating point and in matrix inversion at each frequency, so machine-precision agreement is not realistic.

Every component gets a unit test in both classes from the moment it lands. Validation is not a separate phase — it is built into the development loop. If a component cannot be tested against Finesse, it cannot be merged.

## 8. Phased Plan

### Phase 1 — Directed connector spike and one component
- Define `OpticalOutput` and `OpticalInput` directed sub-port connectors carrying the real and imaginary parts of a complex envelope amplitude.
- Implement `Beamsplitter` (the most general 2-port linear scatterer) with directed ports.
- Verify that `mtkcompile` produces the expected equations and that two beamsplitters can be wired face-to-face cleanly.
- **This is the make-or-break milestone.** If the directed-connector design produces ugly or wrong equations under `mtkcompile`, the architecture must be reconsidered before going further.

### Phase 2 — Plane-wave optomechanics
- Add `Laser`, `Space`, `Mirror`, and `Photodetector` (sensor with zero loading, modeled on `VoltageSensor` in `MTKStdLib/src/Electrical/Analog/sensors.jl:85-105`).
- Build a Fabry-Perot cavity as a reusable `@mtkmodel` subassembly.
- Solve the DC operating point and verify cavity buildup matches the analytical formula.
- Implement `to_finesse` for each component and validate the cavity against Finesse at `rtol=1e-8`.
- **At the end of this phase, every node represents a single complex amplitude (a 1×1 scalar).** This phase covers TEM₀₀ plane-wave optics only.

### Phase 3 — Two-quadrature optomechanics
- Promote each optical node from a single complex amplitude to a 2×2 representation with both quadratures (amplitude and phase, or upper and lower sideband).
- This is what unlocks linearized small-signal optomechanics: radiation pressure couples the amplitude quadrature to mechanical motion, which couples back to the phase quadrature.
- Add radiation-pressure coupling between optical power and a `Translational.MechanicalPort` on `Mirror` and `Beamsplitter`.
- Wire up a suspended cavity using MTK Standard Library masses and springs.
- Validate against Finesse's quadrature-resolved analyses.

### Phase 4 — AC analysis path
- Use `linearize_symbolic` to extract `(A, B, C, D)` for an MTK optical model with a designated input (laser amplitude, laser phase, mirror position) and output (photodetector signal).
- Implement `H(iω) = C(iωI − A)⁻¹B + D` evaluation over a frequency sweep.
- Validate against Finesse's `frequency_response` action at `rtol=1e-3`.
- This is where the hybrid design pays off: a single MTK model definition serves both DC operating point and AC sideband analysis.

### Phase 5 — Closed-loop electro-opto-mechanics
- Add control electronics from MTK Standard Library `Electrical` and `Blocks`.
- Demonstrate a closed-loop length-locked cavity, validating loop gain against a Finesse `Lock` action.
- Demonstrate digital control with sampled-data feedback using MTK's discrete-time machinery, mixed with continuous optical and mechanical dynamics.

## 9. Risks

In rough priority order:

1. **Directed connector behavior under `mtkcompile`.** The Phase 1 spike resolves this. If `mtkcompile` produces sensible equations from face-to-face directed connections, the rest of the plan follows. If not, the architecture must be reconsidered.
2. **Operating-point convergence.** Optical interferometers have many near-degenerate operating points (multiple cavity longitudinal modes). MTK's nonlinear solver may need careful initial guesses. Finesse handles this by solving the carrier first; we don't get that for free.
3. **Complex amplitudes as real pairs.** Workable but doubles equation count and obscures formulas. Worth checking whether `Symbolics.@variables x(t)::Complex` survives through `mtkcompile`; if it does, prefer it.
4. **Round-trip exporter parameter mapping.** Finesse's parameter validation is strict. Every component needs an exporter test from day one.
5. **Performance ceiling at scale.** OpticalModelingToolkit will not compete with Finesse on LIGO-scale models with thousands of HOMs. We compete on composability and electro-mechanical coupling, not on raw matrix size.
6. **Higher-order modes** require a fundamentally different approach (see §10) and should not be attempted under the current architecture.
7. **Quantum noise** is a research project, not a feature.

## 10. Deferred — Higher-Order Modes as a Migration Investigation

Finesse's most distinctive feature is its Hermite-Gaussian mode expansion: every optical node carries an array indexed by `(n, m)` HOM pairs, beam parameters propagate via ABCD matrices, and components have HOM-coupling matrices computed by Cython kernels in `finesse3/src/finesse/components/modal/*.pyx`. This is intrinsically a per-frequency, per-mode-pair linear algebra problem that does not decompose cleanly into MTK symbolic equations.

Studying this system and assessing what a Julia migration would look like is its own interesting investigation, separate from the OpticalModelingToolkit core. Open questions for that investigation:

- Where does the Hermite-Gaussian basis enter the Finesse equations, and is the basis-dependence purely a per-component coupling matrix or does it interact with the network topology?
- Could the mode-coupling matrices be precomputed once per geometry change and then plugged into MTK as `@register_symbolic` opaque functions, leaving the per-mode field amplitudes as MTK variables?
- How is beam tracing (`finesse3/src/finesse/tracing/`) coupled to the equation-building pipeline, and could it run as a non-symbolic pre-pass before MTK sees the model?
- Is there a numerical core that could be reimplemented in Julia (LinearSolve.jl + KLU) without needing the full Finesse symbolic and KatScript machinery on top?

These questions are out of scope for Phase 1–5. They are recorded here so they aren't lost.
