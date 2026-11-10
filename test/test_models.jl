# # Shared models: one file, two simulators
#
# Each file in `models/` is a KatScript model that is the single source of truth
# for both codes. Finesse.jl solves its frequency response. OMT builds the same
# network as a time-domain DAE with [`build_kat`](@ref): optics, suspension and
# servo in one ModelingToolkit model.
#
# This page checks the two against each other without Finesse installed.
# Finesse's responses are stored in `models/reference/`; `test/finesse/reference.jl`
# writes them, and the cross-check suite in `test/finesse/` fails if they go stale.
# Three comparisons:
#
# 1. OMT's linearized frequency response against Finesse's;
# 2. OMT's time-domain response to a sinusoidal force, demodulated, against the
#    same Finesse numbers;
# 3. what only the time domain can show: a mirror running away without its servo,
#    and pulled back in with it.

using Test
using OpticalModelingToolkit
using ModelingToolkit
using ControlSystemsBase: ss, poles
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using SciMLBase

const MODELS = joinpath(pkgdir(OpticalModelingToolkit), "models")
kat(model) = read(joinpath(MODELS, "$model.kat"), String)

"Finesse.jl's stored responses of `model`: frequencies [Hz], and one complex vector per response."
function finesse_reference(model)
    path = joinpath(MODELS, "reference", "$model.txt")
    rows = [parse.(Float64, split(l)) for l in eachline(path) if !startswith(l, '#')]
    n = (length(rows[1]) - 1) ÷ 2
    return [r[1] for r in rows], [[complex(r[2k], r[2k + 1]) for r in rows] for k in 1:n]
end

"OMT's frequency response from the force on the suspended mirror to each of `outs(h)`."
function omt_response(src, fs, outs; rep = PadeDelay(1))
    plant, h = build_kat(src; rep)
    o = outs(h)
    mats, _ = linearize_proper(plant, [h.F_z], o; op = Dict(h.z => 0.0, h.v => 0.0))
    G = frequency_response(mats, fs)
    return [G[k, 1, :] for k in eachindex(o)]
end

relerr(a, b) = abs.(a .- b) ./ abs.(b)

# ## The suspended cavity: an optical spring
#
# `models/suspended_fp.kat` is a cavity whose end mirror hangs on a 20 Hz
# pendulum, tuned half a degree past resonance. There the radiation pressure
# rises as the mirror moves in, so the light is a spring, and 50 W makes it much
# stiffer than the suspension.
#
# The spring's damping comes from the cavity's storage time, so a delay
# representation is needed to get it right. `PadeDelay(1)` matches Finesse's
# exact delays; the adiabatic `Algebraic` misses the damping, by up to
# ``4\times10^{-5}`` near the resonance.
#
# The responses are evaluated with [`frequency_response`](@ref), in extended
# precision: with delay states the model is so stiff that a Float64 evaluation
# loses up to ten digits.

fs, (fin_sus,) = finesse_reference("suspended_fp")
omt_sus = only(omt_response(kat("suspended_fp"), fs, h -> [h.z]))
alg_sus = only(omt_response(kat("suspended_fp"), fs, h -> [h.z]; rep = Algebraic()))
@test maximum(relerr(omt_sus, fin_sus)) < 1.0e-9
@test maximum(relerr(alg_sus, fin_sus)) > 20 * maximum(relerr(omt_sus, fin_sus))

# The resonance, which the bare suspension has at 20 Hz:

f_peak = fs[argmax(abs.(fin_sus))]
@test 50 < f_peak < 70
f_peak

# Every fifth frequency: Finesse's ``|z/F|`` [m/N], and the relative error of OMT
# with `PadeDelay(1)` and with `Algebraic`.

[fs abs.(fin_sus) relerr(omt_sus, fin_sus) relerr(alg_sus, fin_sus)][1:5:end, :]

# ## The locked cavity: optics, suspension and servo
#
# `models/locked_cavity.kat` is tuned to the other side of resonance, where the
# optical spring is negative. Without feedback the mirror is unstable. A readout
# of the transmitted power and a lead filter push it back, with a unity-gain
# frequency near 100 Hz.
#
# Open the loop by deleting the actuator line, and the linearized plant has a pole
# in the right half-plane:

const ACT = "wire act servo.p2.o m2.mech.F_z"
open_loop = replace(kat("locked_cavity"), ACT => "")
ol, ho = build_kat(open_loop; rep = PadeDelay(1))
olmats, _ = linearize_proper(ol, [ho.F_z], [ho.z]; op = Dict(ho.z => 0.0, ho.v => 0.0))
growth = maximum(real, poles(ss(olmats...)))
@test growth > 10
growth

# With the loop closed, the response from a force on the mirror to its
# displacement and to the transmitted power matches Finesse's:

fs, (fin_z, fin_P) = finesse_reference("locked_cavity")
omt_z, omt_P = omt_response(kat("locked_cavity"), fs, h -> [h.z, h.readouts[:trans]])
@test maximum(relerr(omt_z, fin_z)) < 1.0e-9
@test maximum(relerr(omt_P, fin_P)) < 1.0e-9
[fs abs.(fin_z) relerr(omt_z, fin_z) abs.(fin_P) relerr(omt_P, fin_P)][1:5:end, :]

# ## Time domain against frequency domain
#
# Now the same file as a time-domain simulation. The force on the mirror is
# ``A\sin(2\pi f t)``, with the amplitude and frequency as parameters, so one
# compiled system serves every drive. The simulation starts from the steady state
# with the drive off ([`steady_state`](@ref)), so the cavity is already full and
# the loop already locked.

src = kat("locked_cavity")
dp, dh = build_kat(src; drive = (0.0, 1.0), rep = PadeDelay(1))
sys = mtkcompile(dp)
settled = steady_state(ODEProblem(sys, [dh.z => 0.0, dh.v => 0.0], (0.0, 0.3)))

# After the transient (the slowest closed-loop pole decays at about 250/s), the
# displacement is ``z = A\,\mathrm{Im}(H e^{i\omega t})``. Demodulating over whole
# periods gives back ``H``:

function demodulate(sol, x, f, T, n)
    ts = range(T - n / f, T; length = 400n + 1)
    w = sol.(ts; idxs = x) .* cis.(-2π * f .* ts)
    integral = step(ts) * (sum(w) - (w[1] + w[end]) / 2)
    return 2im * integral / (n / f)
end

A = 1.0e-8                                           # [N]; z stays below 1e-12 m, well inside the linear range
idx = [11, 16, 21]                                   # 10, 31.6 and 100 Hz
H_td = map(fs[idx]) do f
    prob = remake(settled; p = [dh.drive.drive_A => A, dh.drive.drive_f => f])
    sol = solve(prob, Rodas5P(); reltol = 1.0e-10, abstol = 1.0e-18, initializealg = BrownFullBasicInit())
    @test SciMLBase.successful_retcode(sol)
    demodulate(sol, dh.z, f, 0.3, max(1, floor(Int, 0.15f))) / A
end
@test maximum(relerr(H_td, fin_z[idx])) < 1.0e-5

# The frequency, Finesse's ``|H|`` and phase [deg], and the same from the
# time-domain simulation:

[fs[idx] abs.(fin_z[idx]) rad2deg.(angle.(fin_z[idx])) abs.(H_td) rad2deg.(angle.(H_td))]

# ## What only the time domain can show
#
# Finesse linearizes about the operating point, and so answers questions about
# small motions around a working lock. Everything else is a time-domain question.
# Here are two: the open-loop mirror nudged by 0.1 pm runs away, and the locked
# mirror, started a fifth of a linewidth off, is pulled back.

linewidth = 1064.0e-9 / (4π) * (1 - sqrt(0.99 * 0.999)) / (0.99 * 0.999)^(1 / 4)

function nudge(src, z0; T)
    p, h = build_kat(src; drive = (0.0, 1.0), rep = PadeDelay(1))
    s = mtkcompile(p)
    prob = steady_state(ODEProblem(s, [h.z => 0.0, h.v => 0.0], (0.0, T)))
    prob = remake(prob; u0 = [h.z => z0])
    sol = solve(prob, Rodas5P(); reltol = 1.0e-9, abstol = 1.0e-20, initializealg = BrownFullBasicInit())
    return sol, h
end

## open loop: the anti-spring wins
sol_ol, h_ol = nudge(open_loop, 1.0e-13; T = 0.1)
@test SciMLBase.successful_retcode(sol_ol)
@test abs(sol_ol[h_ol.z][end]) > 100 * 1.0e-13

## closed loop: the servo pulls the mirror in
z0 = 0.2 * linewidth
sol_cl, h_cl = nudge(src, z0; T = 0.1)
@test SciMLBase.successful_retcode(sol_cl)
@test abs(sol_cl[h_cl.z][end]) < 1.0e-3 * z0

# Where each mirror ends up after 0.1 s, as a fraction of where it started:

(open_loop = sol_ol[h_ol.z][end] / 1.0e-13, locked = sol_cl[h_cl.z][end] / z0)
