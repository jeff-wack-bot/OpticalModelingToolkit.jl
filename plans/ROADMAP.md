# Roadmap

The current state, what's open, and what was decided against. The executed
design plans (`REFINED_PLAN.md`, `PLAN_rationalize_late.md`) live in git
history; `docs/src/design.md` records the decisions that came out of them.

## Done

- Directed optical ports with `link`; components as complex scattering tables
  (Mirror, MovableMirror, Beamsplitter, FreeSpace, Laser, LaserModulated,
  Vacuum, Photodetector).
- Propagation representations `Algebraic`/`PadeDelay(n)` (all-pass).
- `linearize_proper`: proper, DC-exact linearization through ControlSystemsBase.
- `to_scattering`: exact transcendental frequency response; the fields match
  Breadboard's `ScatteringSystem`.
- Radiation pressure: optical springs from linearization, and time-domain
  lock acquisition with an MTK servo.
- Subassemblies (`expose`), flattened by `to_scattering`.

## Open

- **Mechanical port.** `MovableMirror` exposes `x` (input) and `F` (output) as
  signals, and mechanics are written as equations. The alternative is an
  acausal `ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition`
  flange, which plugs into standard-library masses and springs at the cost of
  a heavy dependency. Undecided.
- **`to_scattering` at the model's parameter values.** It currently uses
  constructor values, so `remake`d parameters are invisible to it.
- **Upstream MTK issues** found here, worth filing:
  1. `linearize` could eliminate an invertible index-1 algebraic block itself
     (what `linearize_proper` does).
  2. Tearing divides by parameter-dependent coefficients that can be
     numerically zero, making the compiled system silently singular (worked
     around with `irreducible` incoming fields).
  3. Delayed references to array/connector variables, and to observed
     variables, aren't supported. This blocks exact delays in the time
     domain.

## Decided against

- **Finesse as a runtime or test dependency.** The validation surface is
  analytic references, the two-photon A-matrix identity, and representation
  convergence. A Finesse oracle, if wanted, belongs in a separate environment.
  (Finesse has no aarch64 wheel.)
- **`ExactDelay` in the time domain.** Spiked on MTK v11.31 and found a no-go
  (upstream issue 3): exact delays break every algebraic loop, and the result
  is a delay-algebraic recursion outside DelayDiffEq's DDE class.
- **A sparse/KLU numerical core.** These networks are small and dense.
- **Higher-order spatial modes, nonlinear optics, quantum noise as a
  primitive.** Outside a real-valued plane-wave DAE.
