# Design

A handful of decisions shape this package. Each one closes off some
alternatives, and the rest of the code follows from them.

## 1. Optics is one more ModelingToolkit domain

OMTK doesn't build a simulator. It adds components to MTK's, so a mirror, a
suspension and a servo filter are equations in a single model. Everything
downstream belongs to the SciML ecosystem: `mtkcompile`, `NonlinearProblem`,
`ODEProblem`, linearization, ControlSystems. OMTK writes no solvers, sweeps
or plotting. This is why a time-domain lock acquisition with nonlinear
electronics in the loop, which a frequency-domain field solver like Finesse
can't do, comes for free.

The cost is that everything must be a real-valued DAE. That rules out
higher-order spatial modes and nonlinear optics, and makes quantum noise a
separate layer on top of the linearized model.

## 2. A field is a complex envelope, realified at the port

The physics of every component is linear in the complex envelope amplitude, so
components are written in complex arithmetic. MTK is real, so a port carries
the pair `a = [re, im]`. The translation happens in exactly one place (`field`
reads a connector as a `Complex{Num}`, `emit` writes one), so each complex
coefficient becomes a 2×2 real block automatically.

Realification isn't only plumbing. Linearizing the real-pair model around a
carrier gives the two-photon (quadrature) picture directly. The scalar
sideband picture is the same physics in another basis, related by the fixed
matrix ``A = \tfrac12\begin{bmatrix}1 & 1\\ -i & i\end{bmatrix}``. The test
suite checks this identity exactly (see *The transcendental scattering bridge*).

## 3. Ports are directed pairs, wired with `link`

A physical optical port carries light both ways. MTK's `connect` either
equates variables or sums flows, and neither of those means "my outgoing field
is your incoming field". So a port is two directed connectors, `i::OpticalIn`
and `o::OpticalOut`, and `link(p, q)` emits one `connect` per direction. This
needs no custom connection semantics and mirrors Finesse's node model, where
every port has an input and an output node.

## 4. A component is a scattering table

Each component is one function returning `(out_port, in_port) => amplitude`
entries. That function is evaluated twice:

- **symbolically**, to give the component's MTK equations, and
- **numerically**, to give the `OpticalComponent` metadata (a `Scatterer`,
  `Propagator` or `Terminal`) that [`to_scattering`](@ref) reads.

The time-domain model and the exact frequency-domain lowering can't drift
apart, because they come from the same table. A new component written this way
joins both with no changes elsewhere.

## 5. Lengths are macroscopic; positions are tunings

A `FreeSpace` is a pure delay `L/c`. Its carrier phase `kL` is taken to be a
whole number of wavelengths, and microscopic positions are a separate tuning
`ϕ` on each mirror or beamsplitter, which reflects with
``\pm r e^{\pm 2i\phi}``. This is how Finesse and Breadboard describe
interferometers. A 10 cm cavity is resonant at zero tuning with no
sub-wavelength bookkeeping of `L`. Numerically, no phase is ever computed from
a ~10⁶ rad argument whose roundoff the cavity finesse would amplify. A
`MovableMirror`'s position `x` simply adds `kx` to its tuning.

## 6. The field network is solved forward, never inverted

Incoming fields are marked `irreducible`, so they always stay unknowns of the
compiled system. The only variable MTK's tearing may solve for in a component
equation is the outgoing field, whose coefficient is 1.

Without this, tearing is free to solve `o = c·i` for `i`, dividing by `c`. The
coefficient `c` may be `√T` or the `cos 2ϕ`/`sin 2ϕ` of a tuning (including a
`MovableMirror` position), which are numerically zero for a perfect reflector
or at particular tunings, such as a cavity exactly on resonance. The compiled
system is then silently singular. Before this rule, it produced an optical
spring of 2×10⁹ N/m where the true value is 0. The price of forbidding it is a few
more algebraic unknowns, and these networks are small.

## 7. How a delay is represented is chosen per analysis

Free-space propagation is the only place where optical dynamics enter. *How*
the delay is represented is a property of the analysis, not of the component:
`FreeSpace(rep = ...)` takes `Algebraic()` (no states, the adiabatic limit) or
`PadeDelay(n)` (an order-n all-pass cascade). A model builder threads `rep`
down from the top, so switching representation is a one-argument change.

The all-pass choice matters. Every section has ``|H(i\omega)| = 1``, so
passivity holds at every order and only the phase is approximate. A low-pass
"delay" silently dissipates energy above ``1/\tau``.

## 8. Rationalize as late as possible, and in the frequency domain never

A time-domain integrator needs a finite state vector, so `PadeDelay` is as
late as rationalization can happen there. Frequency-domain questions don't
need it at all: [`to_scattering`](@ref) lowers the composed model to the form

```math
G(s) = M_3\,\bigl(E(-s) - M_1\bigr)^{-1} M_2 + M_4,
\qquad E(s) = \mathrm{diag}(e^{-s\tau_k}),
```

with every delay exact, so FSR and comb physics are included. This is also the
ground truth that `PadeDelay(n)` converges to. The form is the one Breadboard's
`ScatteringSystem` uses, so a model can be handed to Breadboard field for
field.

## 9. Linearization eliminates the algebraic block

Circulating-field loops survive compilation as algebraic unknowns (and decision
6 guarantees they do). `ModelingToolkit.linearize` keeps them as states behind
differentiated constraints, which gives a realization with input-derivative
terms that is singular at ``s = 0``. [`linearize_proper`](@ref) instead
eliminates the index-1 algebraic block statically from the Jacobians MTK
already computes. The states are then exactly the differential unknowns, and
the result is proper and exact at DC.

## 10. Tests are the documentation

Every file in `test/` is a [Literate](https://github.com/fredrikekre/Literate.jl)
script: `runtests.jl` runs it, and `docs/make.jl` renders it as an *Examples*
page. Each test states the physics it checks against, analytically where
possible. Inside a code block, write comments with `##`. A line starting with
`# ` becomes prose and splits the block.
