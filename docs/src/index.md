# OpticalModelingToolkit.jl

OpticalModelingToolkit (OMTK) adds an optical domain to
[ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl): mirrors,
beamsplitters, free space, lasers and photodetectors that compose with
`connect` like any other MTK component. They therefore share one model with
mechanics, electronics and control loops, and solve with the SciML stack.

The scope is plane-wave (TEM₀₀) optics on the slowly-varying complex envelope:
DC operating points, small-signal transfer functions, time-domain transients,
and exact frequency-domain responses.

## A Fabry-Perot cavity

```julia
using OpticalModelingToolkit, ModelingToolkit, NonlinearSolve
using ModelingToolkit: t_nounits as t

@named laser = Laser(A = 1.0)
@named m1 = Mirror(R = 0.99)          # reflective face p1, toward the cavity
@named space = FreeSpace(L = 0.1)
@named m2 = Mirror(R = 0.999)         # resonant at zero tuning ϕ
@named pd = Photodetector()

eqs = [link(laser.p1, m1.p2); link(m1.p1, space.p1); link(space.p2, m2.p1); link(m2.p2, pd.p1)]
@named cavity = System(eqs, t; systems = [laser, m1, space, m2, pd])

sys = mtkcompile(cavity)
sol = solve(NonlinearProblem(sys, Dict(u => 0.0 for u in unknowns(sys))))
sol[pd.P]    # transmitted power
```

From the same model:

- **Transfer functions**: `linearize_proper(cavity, inputs, outputs)` returns a
  proper `(A, B, C, D)`. Hand it to `ControlSystemsBase.ss` for `freqresp`,
  `bodeplot`, `poles` and so on.
- **Exact frequency response**: `to_scattering(cavity)` returns the cavity
  with every free-space delay kept as an exact `e^{-sτ}`.
- **Delay dynamics in the time domain**: rebuild with
  `FreeSpace(rep = PadeDelay(n))`.

## One model, two simulators

A model written as a KatScript file in `models/` is read by both OMT and
[Finesse.jl](https://git.mccullerlab.com/software/finesse.jl): Finesse solves its
frequency response, and [`build_kat`](@ref) builds the same network, including
suspensions, readouts and servo filters, for OMT to simulate in time. See
[Shared models with Finesse.jl](models.md).

## Where to read next

- [Design](design.md): the handful of decisions that shape the package.
- [Solver strategy](solvers.md): which SciML problem types and solvers to use,
  and why.
- [Shared models with Finesse.jl](models.md): the single-source-of-truth
  models, the conventions that map one code onto the other, and how they compare.
- [Simulation projects](projects.md): what this combination is good for.
- The worked examples are the test suite. Every page under *Examples* is a test
  file, rendered.
- [API](api.md).
