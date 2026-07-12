# Solver strategy

How OpticalModelingToolkit uses the SciML ecosystem, in one sentence: **emit
the right problem types and let the ecosystem dispatch; build no numerical
core.** (See `plans/PLAN_rationalize_late.md` §4 for the full rationale.)

## DC operating point

The steady-state field network is a `NonlinearProblem` solved with
NonlinearSolve. The field network is linear given the inputs, so Newton
converges in one step; NonlinearSolve already routes its inner solves through
LinearSolve.

Declare `guesses` — never defaults — on all algebraically-determined
component variables. Defaults become `Initial()` constraints and
overdetermine DAE initialization (this bug was found and fixed during the
rtsmodel phase-3 work; `Photodetector.P` and `MovableMirror.x` carry guesses
for exactly this reason).

There is no sparse/KLU numerical core, deliberately: OMTK systems are dense
and small (dozens of states), so LAPACK via the LinearSolve defaults is
correct, and it is already in the stack transitively. Revisit only if models
ever grow by two orders of magnitude.

## AC analysis

The route is: numeric linearization at the operating point, then everything
downstream lives in ControlSystemsBase:

```julia
mats, ssys = linearize_proper(model, [laser.δA, laser.δQ], [pd.P])

using ControlSystemsBase
sys = ss(mats...)
freqresp(sys, ωs)     # or bodeplot, poles, dcgain, …
```

OMTK writes no frequency-sweep code at all, and you land in the same toolbox
Etalon's `pade` output lives in — one plotting/analysis idiom across both
tools.

### Why `linearize_proper` and not `ModelingToolkit.linearize`

The circulating-field loops of a cavity survive `mtkcompile` as torn
algebraic unknowns (tearing cannot eliminate them), and every all-pass
`PadeDelay` section has direct feedthrough, so the loops touch the inputs.
`ModelingToolkit.linearize` keeps those algebraic unknowns as states by
differentiating their constraints, which

- demands input-derivative columns (`allow_input_derivatives = true`,
  `B = [B_u  B_u̇]`), and
- leaves one spurious zero eigenvalue per algebraic unknown, so the returned
  realization is singular at `s = 0` — `dcgain` is wrong — and
  `ControlSystemsMTK.named_ss`'s descriptor-based cleanup (`dss2ss`) returns
  wrong frequency responses for these systems.

For an index-1 DAE the honest reduction is static elimination of the
algebraic block, `z = -g_z⁻¹(g_x x + g_u u)`, which is what
[`linearize_proper`](@ref) does with the Jacobians MTK already computes. The
result is proper (no `u̇`), regular at `s = 0`, and has exactly the
differential unknowns as states. `ModelingToolkit.linearize_symbolic` remains
useful as a pedagogical oracle (see `test/test_phase1b_fabry_perot_ac.jl`),
with the same caveat about validity at `s = 0`.

### Exact frequency-domain answers

For frequency-domain questions you need not rationalize at all: lower the
model with [`to_scattering`](@ref) and evaluate the exact transcendental
response, FSR/comb physics included. With Etalon loaded,
`Etalon.ScatteringSystem(model)` gives `freqresp` and `pade(sys; order)`
directly.

## Time domain

The composed optical network compiles to a **mass-matrix DAE**: the algebraic
optical constraints survive `mtkcompile`, and tearing does not eliminate the
circulating-field loops. Therefore:

- Use Rosenbrock-family solvers: `Rodas5P` or `Rosenbrock23`. Explicit RK
  refuses mass matrices.
- With per-tick parameter callbacks (RTSCosim co-simulation), pass
  `initializealg = BrownFullBasicInit()`. MTK v11's default re-initialization
  after callbacks re-enforces `Initial(x)` and silently pins DAE plant states
  (the rtsmodel phase-3 trap; RTSCosim's `simulate` already defaults to
  this).
- Stiffness is a non-issue: `PadeDelay` states live at `θ = τ/n` (~ns) and
  callback tstops cap the step size anyway; Rosenbrock methods jump straight
  onto the slow manifold, at no measured cost versus the algebraic
  representation.

## Choosing a propagation representation

The representation is an analysis-time choice threaded in from the top of a
model — never hard-coded in a component:

| rep | states per `FreeSpace` | use for |
|-----|------------------------|---------|
| `Algebraic()` | 0 | audio-band studies of short cavities (adiabatic limit) |
| `PadeDelay(n)` | 2n | anything approaching `c/L`; FSR physics needs `n ≳ ωτ` |
| `to_scattering` | — (exact, frequency domain only) | FSR/comb physics, convergence ground truth |

`PadeDelay` sections are all-pass: passivity holds exactly at every order,
and only the phase carries the `(ωθ)³/12` per-section error, `∝ n⁻²` at fixed
`ωτ`. If a study ever wants `n > ~20`, the honest answer is the exact
scattering bridge, not more states.
