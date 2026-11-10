# Time the shared models in both codes.
#
#     julia --project=test/finesse test/finesse/benchmark.jl
#
# Every number is the second of two runs, so compilation is excluded except where a
# step compiles by design (`mtkcompile` builds and compiles a new model each time).

using Test
using LinearAlgebra
using OpticalModelingToolkit
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using NonlinearSolve
using SciMLBase
using ControlSystemsBase: ss, freqresp
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve

include("common.jl")
include("kat_models.jl")

"Median seconds of `n` runs of `f`, after one warm-up run."
function timed(f; n = 3)
    f()
    return sort([@elapsed(f()) for _ in 1:n])[(n + 1) ÷ 2]
end

dense = exp10.(range(0, 3; length = 1001))
src = kat_source("locked_cavity")

## Finesse: parse, build and solve the closed-loop response
t_fin31 = timed(() -> finesse_response(src, REFERENCE_FREQS, "m2.mech.F_z", "m2.mech.z"))
t_fin1001 = timed(() -> finesse_response(src, dense, "m2.mech.F_z", "m2.mech.z"))

## OMT: build the model and linearize it (compiles), then evaluate the response
function omt_linearize()
    plant, h = build_kat(src; rep = PadeDelay(1))
    mats, _ = linearize_proper(plant, [h.F_z], [h.z]; op = Dict(h.z => 0.0, h.v => 0.0))
    return ss(mats...)
end
t_omt_build = timed(omt_linearize; n = 1)
sys = omt_linearize()
mats = let plant_h = build_kat(src; rep = PadeDelay(1))
    first(linearize_proper(plant_h[1], [plant_h[2].F_z], [plant_h[2].z]; op = Dict(plant_h[2].z => 0.0, plant_h[2].v => 0.0)))
end
t_omt31 = timed(() -> frequency_response(mats, REFERENCE_FREQS))
t_omt1001 = timed(() -> frequency_response(mats, dense))
t_csb1001 = timed(() -> freqresp(sys, 2π .* dense))

## OMT time domain: one 0.3 s simulation of the driven loop, from a compiled system
dp, dh = build_kat(src; drive = (0.0, 1.0), rep = PadeDelay(1))
csys = mtkcompile(dp)
settled = steady_state(ODEProblem(csys, [dh.z => 0.0, dh.v => 0.0], (0.0, 0.3)))
prob = remake(settled; p = [dh.drive.drive_A => 1.0e-8, dh.drive.drive_f => 31.6])
t_td = timed(() -> solve(prob, Rodas5P(); reltol = 1.0e-10, abstol = 1.0e-18, initializealg = BrownFullBasicInit()))

println("| step | time |")
println("|:-----|-----:|")
for (name, s) in (
        ("Finesse.jl: frequency response, 31 points", t_fin31),
        ("Finesse.jl: frequency response, 1001 points", t_fin1001),
        ("OMT: `build_kat` + `linearize_proper` (compiles the model)", t_omt_build),
        ("OMT: `frequency_response` (BigFloat), 31 points", t_omt31),
        ("OMT: `frequency_response` (BigFloat), 1001 points", t_omt1001),
        ("OMT: `ControlSystemsBase.freqresp` (Float64, inexact here), 1001 points", t_csb1001),
        ("OMT: 0.3 s time-domain simulation of the loop", t_td),
    )
    println("| $name | $(s < 1.0e-3 ? string(round(1.0e6 * s; sigdigits = 3), " µs") : s < 1 ? string(round(1.0e3 * s; sigdigits = 3), " ms") : string(round(s; sigdigits = 3), " s")) |")
end
