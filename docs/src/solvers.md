# Solver strategy

OMTK emits standard SciML problem types and lets the ecosystem dispatch. It
has no numerical core of its own.

## DC operating point

Compile, then solve a `NonlinearProblem`:

```julia
sys = mtkcompile(model)
sol = solve(NonlinearProblem(sys, Dict(u => 0.0 for u in unknowns(sys))))
```

The unknowns are the incoming fields of every port (they are `irreducible`,
see [Design](design.md) §5). Given its inputs, the field network is linear,
so zero is a fine starting point and Newton converges in one step. A
`MovableMirror` position makes it nonlinear in `x`, but `x` is bound by an
equation, so the field solve stays linear.

Component variables that are algebraically determined, such as
`Photodetector.P` and `MovableMirror.F`, carry *guesses*, never defaults. A
default becomes an `Initial()` constraint and overdetermines initialization.
Follow the same rule in your own components.

There is no sparse or KLU path, by design. These systems are small and dense,
so the LinearSolve defaults are right.

## Transfer functions

Linearize at an operating point, then work in ControlSystemsBase:

```julia
mats, _ = linearize_proper(model, [laser.δA, laser.δQ], [pd.P]; op = Dict(...))

using ControlSystemsBase
sys = ss(mats...)
freqresp(sys, ωs)     # or bodeplot, poles, dcgain, …
```

Inputs are `[input = true]` variables, such as `LaserModulated.δA`/`δQ` or
`MovableMirror.x`. Pass the value of a nonlinear input such as `x` in `op`,
since the linearization depends on it.

Use [`linearize_proper`](@ref), not `ModelingToolkit.linearize`. The cavity's
circulating fields are algebraic unknowns, and every `PadeDelay` section has
direct feedthrough, so the field loops touch the inputs.
`ModelingToolkit.linearize` keeps those unknowns as states by differentiating
their constraints. This

- demands input-derivative columns (`allow_input_derivatives = true`), and
- leaves one spurious zero eigenvalue per algebraic unknown, so the
  realization is singular at `s = 0` and `dcgain` is wrong.
  `ControlSystemsMTK.named_ss`'s descriptor cleanup returns wrong frequency
  responses for these systems.

`linearize_proper` eliminates the index-1 algebraic block statically,
``z = -g_z^{-1}(g_x x + g_u u)``, using the Jacobians MTK already computes. The
result is proper, regular at `s = 0`, and its states are exactly the
differential unknowns.

## Exact frequency response

For frequency-domain questions, skip rationalization entirely and use
[`to_scattering`](@ref):

```julia
parts = to_scattering(model)
G(s) = parts.M3 * ((Diagonal(exp.(s .* parts.delays)) - parts.M1) \ parts.M2) + parts.M4
```

The channels `parts.u` and `parts.y` are named after the terminals (lasers,
vacua, photodetectors). With Etalon loaded, `Etalon.ScatteringSystem(model)`
wraps the same data.

## Time domain

The compiled model is a mass-matrix DAE: the field network stays algebraic.
Therefore:

- Use Rosenbrock-family solvers, `Rodas5P` or `Rosenbrock23`. Explicit RK
  methods refuse mass matrices.
- If parameter-changing callbacks run during integration (for example
  co-simulation with a digital controller), pass
  `initializealg = BrownFullBasicInit()`. MTK's default re-initialization
  after a callback re-enforces `Initial(x)` and silently pins the plant
  states.
- `PadeDelay` states live on the `τ/n` timescale (nanoseconds). Stiffness is
  not a problem for Rosenbrock methods, which step straight onto the slow
  manifold.

## Choosing a propagation representation

| rep | real states per `FreeSpace` | use for |
|-----|-----------------------------|---------|
| `Algebraic()` | 0 | audio-band studies of short cavities (adiabatic limit) |
| `PadeDelay(n)` | 4n | anything approaching `c/L`; FSR physics needs `n ≳ ωτ` |
| `to_scattering` | (exact, frequency domain only) | FSR/comb physics, ground truth for convergence |

The phase error of `PadeDelay(n)` is ``(\omega\tau)^3/12n^2``. If a study needs
`n > ~20`, use the exact scattering form rather than more states.
