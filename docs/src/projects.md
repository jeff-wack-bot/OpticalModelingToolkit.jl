# Simulation projects

OMT's niche is the time domain of hybrid systems: light, suspensions and
electronics in one model, with the nonlinearities left in. Paired with Finesse.jl
on a shared model, a linear design can be checked against what actually happens
when the loop is closed, opened, saturated or started from the wrong place. These
are projects that use exactly that, roughly in order of how much they need
beyond what OMT has now.

## Within reach now

**Lock acquisition of an arm cavity.** A LIGO arm is a 4 km cavity of finesse
~450 whose mirrors swing by microns while the linewidth is nanometres; the light
flashes as a mirror crosses resonance, and a servo has to catch it. The questions
are about capture range and the guided-lock tricks that extend it, and they are
nonlinear and transient by nature. `locked_cavity.kat` already acquires from a
fifth of a linewidth; the project is a realistic sweep velocity, a realistic
servo, and the capture probability as a function of both.

**Resonance sweeps and cavity ringing.** A cavity swept through resonance faster
than its storage time rings: the transmitted power oscillates after the peak,
with a chirp set by the sweep rate and the delay. It is a classic measurement of
the storage time and the finesse, it has an analytic answer, and it tests the
delay representation directly, since `Algebraic` cannot ring and `PadeDelay(n)`
must converge on the exact result.

**Actuator saturation and lock loss.** Real drivers saturate and DACs have a
range. Put a saturation in the servo (an MTK block, which Finesse cannot
represent), drive the plant with a disturbance that the linear design rejects,
and find where the loop breaks and whether it recovers.

**The optical-spring instability, stabilized and unstabilized.** The locked
cavity's anti-spring is unstable without feedback; with feedback, how much gain
margin there is depends on power. Sweep the input power and map where the loop
fails, in both codes: the closed-loop poles from Finesse's response and OMT's
linearization, and the onset of runaway from OMT's time domain.

## A feature or two away

**Sampled-data control, as LIGO runs it.** LIGO's servos are digital, sampled at
16 or 65 kHz with ADC and DAC delays and anti-aliasing filters. ModelingToolkit
has clocked, discrete-time variables; a ZPK becomes a biquad cascade on a clock,
and the model becomes a hybrid continuous-discrete system. The project is a
loop whose design was done in the continuous domain, simulated as it actually
runs, with the phase lost to sampling and delay showing up as reduced margin.

**Multi-stage suspensions and hierarchical actuation.** A LIGO test mass hangs
from a quadruple pendulum, actuated at several stages with crossovers between
them, and damped by local sensors at the top. In OMT a suspension is equations,
so a four-stage chain with its local damping loops and the length loop on top is
one model; the KatScript side needs `suspension_zpk` (or a reduced model of the
chain) for Finesse to see the same plant.

**RF sidebands and Pound-Drever-Hall locking.** Most interferometer locks are
PDH: phase-modulate the laser at RF, demodulate the reflected light. OMT's fields
are a single carrier envelope; adding a few sideband envelopes per port (carrier
and ``\pm\Omega``) would make PDH, and the sensing schemes built on it,
expressible, and the shared models could then include Finesse's `modulator` and
`readout_rf`.

## Further out

**A coupled-cavity interferometer.** Power and signal recycling make the
interferometer a set of coupled cavities whose lock acquisition is a sequence,
with the arms held off resonance while the vertex locks. A plane-wave model of
the full interferometer, with a lock sequence written as logic in the model, is
a large but direct extension of the locked cavity.

**Noise in the time domain.** Finesse computes noise spectra from the
linearization. Driving OMT with stochastic forces and laser noise (a stochastic
DAE) and computing the spectra from the simulation checks those budgets where
the linearization stops being valid, near lock loss or with saturation.
