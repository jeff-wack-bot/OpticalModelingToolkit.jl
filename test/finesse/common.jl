# Shared helpers: run a KatScript model in Finesse.jl, and translate OMT
# parameters into Finesse ones.
#
# ## Conventions, and how they map
#
# Finesse is run with `imaginary_transmission=false` throughout, so both codes
# use real transmission. Measured against Finesse.jl (see `test_conventions.jl`):
#
#   - Finesse reflects with `-r e^{2iφ}` on `p1` and `+r e^{-2iφ}` on `p2`; OMT
#     with `+r e^{2iϕ}` on `p1` and `-r e^{-2iϕ}` on `p2`. So **an OMT mirror is
#     a Finesse mirror turned around**: OMT `p1`/`p2` are Finesse `p2`/`p1`, and
#     OMT tuning `ϕ` [rad] is Finesse `phi = -deg(ϕ)`. The beamsplitter is the
#     same, with OMT `p1, p2, p3, p4` at Finesse `p3, p4, p1, p2` (`alpha=0`).
#     [`mport`](@ref) and [`bsport`](@ref) do the translation.
#   - The other way to match the signs is a quarter-wave offset,
#     `phi = deg(ϕ) + 90`, which keeps the port names. It is right at DC and
#     wrong for sidebands: Finesse treats a tuning as a displacement, whose phase
#     scales with optical frequency as `2φ(1 + f/f₀)`, so a 90° offset adds a
#     real sideband phase. OMT applies the same tuning phase to every sideband.
#     Turning the mirror around leaves zero tuning at zero, so the two codes
#     differ only by that scaling of the tunings actually present.
#   - Amplitude detectors (`ad`) read √W, the same as OMT field amplitudes, and
#     a laser `P` emits amplitude `√P` with its `phase` in degrees.
#   - A `space` puts no carrier phase on the field (macroscopic lengths are whole
#     wavelengths), the same as OMT's `FreeSpace`, so only the sidebands see `L`.
#   - The upper signal sideband at `fsig.f = f` propagates like OMT's
#     `G(s)` at `s = 2πi f`, and `pd1 ... fsig` demodulates with the same phase
#     convention as ControlSystemsBase's `freqresp`.

using Finesse: Finesse, compile_script

"Run a KatScript model in Finesse.jl and return the `Solution`."
finesse_run(src) = Finesse.run(compile_script(src).model)

"The complex single-point reading of detector `name`."
reading(sol, name) = sol[name][1]

"Finesse tuning [deg] equivalent to OMT tuning `ϕ` [rad] on a turned-around surface."
kat_phi(ϕ) = -rad2deg(ϕ)

"Finesse port of OMT mirror port `k`, as `\"name.pN\"`."
mport(name, k) = "$name.p$((2, 1)[k])"

"Finesse port of OMT beamsplitter port `k`, as `\"name.pN\"`."
bsport(name, k) = "$name.p$((3, 4, 1, 2)[k])"

"KatScript for an OMT `Mirror(R, T, ϕ)`; address its ports with [`mport`](@ref)."
kat_mirror(name; R, T = 1 - R, ϕ = 0.0) =
    "mirror $name R=$R T=$T L=$(max(0.0, 1 - R - T)) phi=$(kat_phi(ϕ)) imaginary_transmission=false"

"KatScript for an OMT `Beamsplitter(R, T, ϕ)`; address its ports with [`bsport`](@ref)."
kat_bs(name; R, T = 1 - R, ϕ = 0.0) =
    "bs $name R=$R T=$T L=$(max(0.0, 1 - R - T)) phi=$(kat_phi(ϕ)) alpha=0 imaginary_transmission=false"

"Finesse's default carrier frequency `f₀ = c/λ₀` with `λ₀ = 1064 nm`."
const F0 = Finesse.C_LIGHT / 1064.0e-9

"""
    finesse_freqresp3(src, f, input, outputs) -> Vector{ComplexF64}

Upper-sideband optical transfer functions at frequency `f` from the field
entering `input` (a node such as `"m1.p2.i"`) to each output node.
"""
function finesse_freqresp3(src, f, input, outputs)
    outs = join(("[$o, fsig.f]" for o in outputs), ", ")
    sol = finesse_run(src * "\nfsig(1)\nfreqresp3($f, [[$input, fsig.f]], [$outs])\n")
    return vec(sol.data)
end

"Complex field on an OMT `OpticalIn`/`OpticalOut` in a solved system."
cfield(sol, port) = complex(sol[port.a[1]], sol[port.a[2]])

"Solve an algebraic OMT model's DC operating point from zero."
function omt_dc(model)
    sys = mtkcompile(model)
    sol = solve(NonlinearProblem(sys, Dict(u => 0.0 for u in unknowns(sys))), NewtonRaphson())
    @test SciMLBase.successful_retcode(sol)
    return sol
end

"""
    omt_tf(parts, f) -> Matrix{ComplexF64}

Exact transfer matrix of a `to_scattering` result at frequency `f` [Hz].
"""
function omt_tf(parts, f)
    s = 2π * im * f
    return parts.M3 * ((Diagonal(exp.(s .* parts.delays)) - parts.M1) \ parts.M2) + parts.M4
end
