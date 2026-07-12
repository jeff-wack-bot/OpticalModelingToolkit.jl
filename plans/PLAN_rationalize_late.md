# Plan: rationalize late & the SciML solver strategy

**Status: executed 2026-07-10 — see §10 for outcomes and deviations.**

Response to the findings of Breadboard's port-shape study
(`~/dev/Breadboard/plans/PLAN_port_shapes.md`, Part 5), which examined OMTK
alongside sflu-optomechanics, Finesse 3, and Breadboard/Etalon. Two questions
drive this plan:

1. **How can OMTK rationalize late, in the spirit of Etalon?**
2. **What is the best way for OMTK to use the SciML ecosystem's solvers?**

## 1. The findings this plan answers

From the study (§5.3), plus one defect found while planning:

- **OMTK rationalizes too early.** `LengthDelay` bakes a first-order
  rational approximation of free-space propagation into each component at
  *definition* time, before the network is composed. Cavity FSR physics is
  unrepresentable and the approximation order cannot be swept. Etalon's
  architecture is the counterpoint: the composed `ScatteringSystem` keeps
  every delay as an exact transcendental `e^{-sτ}`, and `pade(sys; order)` is
  a *transformation* applied after composition, with order a per-analysis
  choice.
- **Worse, the current `LengthDelay` dissipates.** `τ·ȯ + o = e^{ikL}·i` is a
  *low-pass*: `|H(iω)| = 1/√(1+ω²τ²) < 1`. True propagation is *all-pass*.
  Every `LengthDelay` instance silently breaks passivity/energy conservation
  above `1/τ` — an error in the physics, independent of the
  rationalize-late question, and fixable at any order (§3).
- **OMTK's real-pair envelope is the realification functor in the wild**
  (study §5.3 obs. 1): linearizing the real-pair DAE around a solved carrier
  already produces two-photon-structured blocks. This makes a specific
  cross-tool validation possible (§5).
- The hand-rolled AC loop (`linearize_symbolic` → `fixpoint_sub` →
  `C(iωI−A)⁻¹B+D` per frequency) duplicates machinery the ecosystem provides;
  REFINED_PLAN §10's deferred "LinearSolve.jl + KLU numerical core?" question
  deserves a definite answer (§4).

## 2. The design principle

Adopt the study's P2 for *representations*, not just port shapes: **how
propagation is represented is an analysis-time choice, not an element
identity.** A space between two mirrors is one component; whether it is an
algebraic phase rotation, an order-n rational delay, or an exact delay is the
analysis's resolution choice — threaded in from the top of a model, exactly
as `mlib`/`shape=` is in sflu/Breadboard.

Concretely: one `FreeSpace` component replacing the `Length`/`LengthDelay`
pair, parameterized by a representation object:

```julia
abstract type PropagationRep end
struct Algebraic  <: PropagationRep end          # current Length: exact at DC, no dynamics
struct PadeDelay  <: PropagationRep; order::Int; end   # order-n all-pass rational
struct ExactDelay <: PropagationRep end          # true delayed argument (DDE) — §6, gated

@component function FreeSpace(; name, L, λ = 1064e-9, rep::PropagationRep = Algebraic())
```

Model constructor functions accept `rep` and pass it down, so switching a
whole model from adiabatic to order-4 to (eventually) exact is a one-argument
change at the top — "rationalize late" in the sense that matters inside a
time-domain DAE framework: **the approximation is chosen per analysis, at
model-build time, sweepable, and never hard-coded in a component
definition.** (MTK cannot defer it later than problem construction — the
integrator needs a finite state vector; deferral past composition belongs to
the frequency domain, where §5 gets it exactly.)

## 3. Work item A — `FreeSpace` and the all-pass Padé ladder

**A1. `Algebraic()`** — the current `Length` equations, unchanged. This is
the adiabatic limit and stays the right default for audio-band studies of
short cavities (validated at 0.1 m scale in rtsmodel's phase-3 lock demo).

**A2. `PadeDelay(n)`** — the workhorse. Realize `e^{-sτ}` as a cascade of
`n` first-order **all-pass** sections, each approximating `e^{-sθ}`,
`θ = τ/n`:

```
H₁(s) = (1 - sθ/2)/(1 + sθ/2)        state form:  ẋ = -(2/θ)(x - u),  y = 2x - u
```

applied identically to both quadratures, with the static rotation `R(2πL/λ)`
kept as a separate algebraic factor (mirroring Etalon folding rotations into
the `M` matrices while `E(s)` stays a pure scalar delay). Properties that
make this the right realization:

- `|H₁(iω)| ≡ 1` exactly — **passivity holds at every order**, fixing the
  `LengthDelay` defect even at n = 1. Only phase is approximate (error
  `~(ωθ)³/12` per section).
- A cascade of first-order sections is numerically well-conditioned where a
  monolithic order-n Padé companion form is not; it is also literally a
  spatial discretization of the line (n segments of length L/n).
- 2n real states per quadrature pair per direction; FSR physics becomes
  representable once `n ≳ ωτ` across the band of interest — and now that's a
  *convergence knob*, not a wall.

**A3. Migration.** `LengthDelay` is deprecated with a docstring stating the
dissipation defect and pointing at `FreeSpace(rep = PadeDelay(1))` (all-pass,
same state count, correct magnitude). The phase-1b AC test gets a
representation-parameterized analytic reference (§7). 

JEFF - this library has exactly 0 users. always put the code in the simplest final state and clean up all refuse and remains of old code paths.

## 4. Work item B — the SciML solver strategy

The answer to "how should OMTK use SciML solvers" is: **emit the right
problem types and let the ecosystem dispatch; build no numerical core.**
Point by point:

- **DC operating point:** keep `NonlinearProblem` + NonlinearSolve. The
  field network is linear given the inputs, so Newton converges in one step;
  NonlinearSolve already routes its inner solves through LinearSolve.
  Declare `guesses` (never defaults) on all algebraically-determined
  component variables — defaults become `Initial()` constraints and
  overdetermine DAE initialization (bug found and fixed during rtsmodel
  phase 3).
- **The REFINED_PLAN §10 "LinearSolve + KLU core" question — answered no.**
  KLU/sparse machinery is what Finesse needs at 10⁴ unknowns with per-point
  refill (study §4.3-6). OMTK systems are dense and small (dozens of
  states); LAPACK via LinearSolve defaults is correct, and it is already in
  the stack transitively. Revisit only if models ever grow two orders of
  magnitude.
- **AC analysis:** replace the hand-rolled loop with the ecosystem path:
  1. solve DC (`NonlinearProblem`),
  2. **numeric** `ModelingToolkit.linearize` at that operating point
     (keep `linearize_symbolic` for pedagogy and as a test oracle only —
     numeric linearization avoids the `fixpoint_sub` parameter-binding pass
     entirely),
  3. `ControlSystemsMTK.named_ss` → a `ControlSystemsBase.StateSpace`
     (verified: ControlSystemsMTK v2.6.0 co-resolves with MTK v11.31),
  4. `freqresp` / `bodeplot` / `poles` / `dcgain` from ControlSystemsBase —
     OMTK writes no frequency-sweep code at all, and users land in the same
     toolbox Etalon's `pade` output lives in (one plotting/analysis idiom
     across both tools).
- **Time domain:** mass-matrix DAEs (the algebraic optical constraints
  survive `mtkcompile`; tearing does not eliminate circulating-field loops —
  measured in rtsmodel phase 3) ⇒ Rosenbrock-family solvers
  (`Rodas5P`/`Rosenbrock23`); explicit RK refuses mass matrices. With
  per-tick parameter callbacks (RTSCosim co-simulation), pass
  `initializealg = BrownFullBasicInit()` — MTK v11's default re-initialization
  after callbacks re-enforces `Initial(x)` and silently pins DAE plant states
  (the rtsmodel phase-3 trap; RTSCosim's `simulate` already defaults to
  this). Stiffness is a non-issue: with `PadeDelay` states at `θ = τ/n ≈ ns`
  and callback tstops capping steps, Rosenbrock methods jump straight onto
  the slow manifold (measured: no cost vs the algebraic representation).
- **Exact delays (if §6 lands):** `DDEProblem` + `MethodOfSteps` from
  DelayDiffEq.

## 5. Work item C — the transcendental bridge (rationalizing *latest*)

For frequency-domain questions OMTK need not rationalize at all: lower the
optical network to Etalon's representation and inherit its exactness.

- **C1 — `to_scattering(model)`.** Walk the component graph: mirrors and
  rotations contribute static entries to `M1..M4`; each `FreeSpace`
  contributes one exact delay `τ = L/c` per direction (rotation folded into
  the adjacent matrix entries, exactly Etalon's `_collapse_linearspaces`).
  Output an Etalon `ScatteringSystem`: exact transcendental `freqresp`,
  FSR/comb physics for free, and `pade(sys; order)` as the late
  transformation into ControlSystems. Scalar (complex-amplitude) fields
  first; OMTK's quadrature pairs are its realification.
- **C2 — the A-matrix keystone test** (study §2.6-3, §4.3-4). For a passive
  network, OMTK's *linearized two-photon* response must satisfy

  ```
  G_omtk(Ω) == A · diag( G_etalon(i(ω₀+Ω)), conj(G_etalon(i(ω₀−Ω))) ) · A⁻¹
  ```

  blockwise per port pair. This single test cross-validates OMTK's
  realification, its linearization path, `to_scattering`, and Etalon's
  solver against each other — a ground-truth harness independent of hand
  derivations, complementing the existing Finesse oracle. It is *the*
  statement both the sideband and quadrature communities agree on.

Scope note: C is not needed for OMTK's core mission (time-domain transients,
digital control, co-simulation). Its value is (a) the validation web and
(b) making "what order do I need?" measurable against exact truth — which is
precisely the study Etalon's `compare_figure` is built around, now available
to OMTK models.

## 6. Work item D — `ExactDelay` in the time domain (gated prototype)

MTK supports DDEs (`DDEProblem` exists; solved via MethodOfSteps), **but**
delayed references require callable variables (`@variables x(..)`) — probed
on MTK v11.31: a plain connector variable errors with "Sym x(t) is not
callable". OMTK's directed connectors are ordinary array variables, so
`out.a ~ R(φ)·in.a(t-τ)` does not compose through `connect` as-is. Options,
in preference order: (i) a post-flattening rewrite that replaces a marker
component's equations with delayed-argument equations on the flattened
system; (ii) frequency-domain exactness only, via §5 (time domain uses
`PadeDelay`); (iii) an upstream MTK issue about delayed connector variables.
Prototype (i) in a spike before committing; do not block Phases A–C on it.

## 7. Validation ladder

1. **All-pass unit test:** `|H(iω)| = 1` to machine precision for every
   section count, plus phase-error scaling `∝ n⁻²` at fixed ωτ (cascade of n
   sections, each error (ωτ/n)³ → total ∝ n⁻²).
2. **FP convergence:** transmission of a Fabry-Perot with `PadeDelay(n)`
   against the exact analytic response `t₁t₂e^{iΩτ}/(1 + r₁r₂e^{2i(kL+Ωτ)})`:
   error decreases monotonically with n over a fixed band; FSR-peak location
   correct to tolerance for n above the rule-of-thumb threshold. (The
   existing phase-1b analytic reference `(1+iωτ)²`-form is the n=…-specific
   one; it becomes representation-parameterized.)
3. **Representation invariance in-band:** `Algebraic` vs `PadeDelay(n)` vs
   §5 exact agree in the audio band for a 0.1 m cavity.
4. **AC-path regression:** ControlSystemsMTK route reproduces the current
   hand-rolled loop's numbers on the phase-1b system (then the loop is
   demoted to an oracle).
5. **A-matrix cross-tool test** (§5 C2) on an undetuned and a detuned FP.
6. **Co-sim smoke:** the rtsmodel cavity-lock demo under `PadeDelay(2)` —
   lock acquisition and measured loop TF unchanged (audio band ≪ 1/τ).

## 8. Phasing

- **Phase A** — `FreeSpace` + `Algebraic`/`PadeDelay`, all-pass sections,
  `LengthDelay` deprecation, tests 1–3. Small and self-contained; fixes the
  passivity defect immediately.
- **Phase B** — AC path via numeric `linearize` + ControlSystemsMTK
  (test 4); documentation of the solver guidance in §4 (guesses, Rosenbrock,
  BrownFullBasicInit) in the package docs.
- **Phase C** — `to_scattering` + the A-matrix harness (tests 5). Etalon is
  unregistered: use a `Pkg.develop` path dependency in the test environment,
  or vendor a minimal `ScatteringSystem` (~100 LOC) if that's too awkward.
- **Phase D** — `ExactDelay` spike (§6), go/no-go, then either land it or
  record the negative result and the upstream issue link.

## 9. Risks

- **MTK DDE ↔ connector composability** — verified real (§6); phase-gated.
- **ControlSystemsMTK API churn** under fast-moving MTK v11 — verified
  co-resolution today; pin in Project.toml and smoke-test first in Phase B.
- **High-order conditioning** — the cascade realization avoids companion-form
  blowup; if a study ever wants n > ~20 the honest answer is §5's exact
  bridge, not more states.
- **Etalon coupling** — keep the C1 interface one-directional (OMTK →
  ScatteringSystem struct) so neither package depends on the other's
  internals; the A-matrix test lives in OMTK's test suite with Etalon as a
  test-only dependency.

## 10. Execution record (2026-07-10)

Phases A–C landed; Phase D spiked to a **no-go**. Per review feedback,
`Length`/`LengthDelay` were deleted outright (no deprecation shims).

- **Phase A** — `FreeSpace(L, λ, c, rep)` with `Algebraic()`/`PadeDelay(n)`
  in `src/components.jl`; all-pass cascade as specified in §3. Validation
  1–3 in `test/test_propagation.jl`: |H| = 1 to 1e-12 at every order, phase
  error ratios match n⁻², FP convergence monotone (max rel err 5.6e-2 →
  7e-4 for n = 1 → 16 over 0.6 FSR), FSR peak within 0.3% at n = 16,
  audio-band representation invariance as predicted (algebraic deviates by
  the cavity-pole correction ~1% at 10 kHz for the 0.1 m cavity).
- **Phase B, one deviation** — the `ControlSystemsMTK.named_ss` route is
  numerically broken for these networks: torn circulating-field loops
  survive `mtkcompile` as algebraic unknowns, MTK's `linearize` keeps them
  as states behind differentiated constraints (input-derivative columns +
  one spurious zero eigenvalue each, so the realization is singular at
  s = 0), and the `dss2ss` cleanup returns wrong frequency responses
  (0.65 rel err at the FSR). The fix is static elimination of the index-1
  algebraic block, implemented as `linearize_proper` (`src/analysis.jl`,
  ~25 lines over MTK's own `linearization_function` Jacobians — exact at
  all frequencies including DC, verified to 1e-10 against
  representation-matched analytic references). Everything downstream is
  ControlSystemsBase (`ss`/`freqresp`/`poles`/`dcgain`) as planned;
  ControlSystemsMTK is not a dependency after all. The old hand-rolled loop
  survives only as the oracle in `test_phase1b_fabry_perot_ac.jl`
  (validation 4). Solver guidance documented in `docs/src/solvers.md`.
  Upstream-issue-worthy: MTK could do this reduction itself when g_z is
  invertible.
- **Phase C** — `to_scattering` (`src/scattering.jl`) walks
  `OpticalComponent` metadata that every component now stamps; FreeSpace
  chains collapse into exact delays with rotations folded into consuming
  matrix entries (Etalon's `_collapse_linearspaces`, §5). Returns the
  `ScatteringSystem` fields as a NamedTuple; a package extension
  (`ext/OpticalModelingToolkitEtalonExt.jl`, weakdep on Etalon) provides
  `Etalon.ScatteringSystem(model)` one-directionally. Validation 5 in
  `test/test_scattering.jl`: exact transcendental FP response over two FSRs
  (1e-6, limited by finesse-amplified large-angle trig roundoff of
  φ ≈ 6e5 rad, not physics); the A-matrix keystone holds blockwise on
  undetuned and detuned cavities — verified *exactly* (all Ω) against the
  rep-matched scalar network for both `Algebraic` and `PadeDelay` (the
  identity is exact for any real-rational representation), plus a
  PadeDelay(16) ↔ exact-delay convergence leg; Etalon's solver matches
  `to_scattering` to 1e-12 through the extension when Etalon is dev'd in.
- **Phase D — no-go on MTK v11.31**, negative result recorded per §8:
  1. connector (array) variables cannot appear delayed — `isdelay` only
     recognizes scalar calls `x(t − τ)`, so `a(t − τ)[i]` is invisible
     (confirms the §6 probe, one level deeper);
  2. the post-flattening rewrite (option i) *does* compose using scalar
     alias variables `uin(t) ~ a[1]` + `uin(t − τ)`, and `mtkcompile`
     accepts it (`is_dde = true`);
  3. but exact delays break every algebraic loop, so tearing eliminates
     *all* field unknowns — the cavity becomes a 0-state system whose
     fields are observed functions of history, and MTK's observed codegen
     cannot evaluate delayed observed references (`UndefVarError` for the
     delayed alias inside the generated function). The physics is a pure
     delay-algebraic recursion, outside DelayDiffEq's DDE class.
  Upstream issue to file: (a) delayed array/connector variables, (b)
  delayed references to observed variables in observed codegen. Time
  domain stays on `PadeDelay`; exact delays live in the frequency domain
  via §5 (option ii), which landed. No `ExactDelay` struct ships — it
  would be dead API.
- **Validation 6** — rtsmodel cavity-lock demo updated for `FreeSpace` and
  run under `PadeDelay(2)`: see rtsmodel `examples/cavity_lock.jl`.
