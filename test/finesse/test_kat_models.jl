# # The shared models: OMT against Finesse.jl
#
# Each file in `models/` is read by both codes. Finesse.jl solves its frequency
# response; OMT linearizes the network it builds from the same text. They must
# agree, and the responses stored in `models/reference/` must be what Finesse
# computes now.

include("kat_models.jl")

"OMT's response of `src` from `F_z` to each output `outs(h)`, at `fs` [Hz]."
function omt_response(src, fs, outs; rep = PadeDelay(1))
    plant, h = build_kat(src; rep)
    o = outs(h)
    mats, _ = linearize_proper(plant, [h.F_z], o; op = Dict(h.z => 0.0, h.v => 0.0))
    G = frequency_response(mats, fs)
    return [G[k, 1, :] for k in eachindex(o)]
end

@testset "stored references are current" begin
    for (model, pairs) in REFERENCES
        fs, stored = read_reference(model)
        @test fs ≈ REFERENCE_FREQS rtol = 1.0e-15
        @test stored ≈ finesse_references(model, pairs) rtol = 1.0e-12
    end
end

# ## Suspended cavity: the optical spring
#
# Detuned, the light acts as a spring. Its damping comes from the cavity's
# storage time, which the adiabatic `Algebraic` representation leaves out, so
# only a delay representation reproduces Finesse's exact-delay response. The
# tuning is varied by rewriting the text, so both codes always read the same
# model.
#
# Any `PadeDelay(n)` is exact here: the cavity's one-way delay is 0.33 ns, so the
# phase error of even one section at 100 Hz is ``(\omega\tau)^3/12 \sim 10^{-22}``.
# The model is very stiff, so the responses are evaluated in extended precision
# with `frequency_response`; in Float64 they lose up to ten digits.

const SUS_LINE = "mirror m2 R=0.999 T=0.001 L=0 phi=0.5"

function sus_kat(phi)
    src = kat_source("suspended_fp")
    occursin(SUS_LINE, src) || error("models/suspended_fp.kat no longer has the m2 line")
    return replace(src, SUS_LINE => "mirror m2 R=0.999 T=0.001 L=0 phi=$phi")
end

relerr(a, b) = maximum(abs.(a .- b) ./ abs.(b))

@testset "suspended cavity: optical spring" begin
    fs = [0.5, 2.0, 5.0, 10.0, 20.0, 30.0, 60.0, 100.0]
    for phi in (0.0, 0.5, -0.5)
        src = sus_kat(phi)
        fin = finesse_response(src, fs, "m2.mech.F_z", "m2.mech.z")
        pade = only(omt_response(src, fs, h -> [h.z]))
        @testset "phi = $phi deg" begin
            @test relerr(pade, fin) < 1.0e-9
            ## every order agrees, since every order is exact at these frequencies
            @test relerr(only(omt_response(src, fs, h -> [h.z]; rep = PadeDelay(3))), fin) < 1.0e-9
            if phi != 0
                ## the adiabatic model misses the optical damping
                alg = only(omt_response(src, fs, h -> [h.z]; rep = Algebraic()))
                @test relerr(alg, fin) > 20 * relerr(pade, fin)
            end
        end
    end
end

# ## Locked cavity: optics, suspension and servo
#
# The closed loop, from a force on the mirror to its displacement and to the
# transmitted power, through the servo the file describes.

@testset "locked cavity: closed-loop response" begin
    src = kat_source("locked_cavity")
    fs = REFERENCE_FREQS
    omt = omt_response(src, fs, h -> [h.z, h.readouts[:trans]])
    @test relerr(omt[1], finesse_response(src, fs, "m2.mech.F_z", "m2.mech.z")) < 1.0e-9
    @test relerr(omt[2], finesse_response(src, fs, "m2.mech.F_z", "trans.DC.o")) < 1.0e-9
    ## and the servo's own output, which a Float64 evaluation got 0.8 % wrong
    ol = replace(src, "wire act servo.p2.o m2.mech.F_z" => "")
    y = only(omt_response(ol, fs, h -> [h.filters[:servo].y]))
    @test relerr(y, finesse_response(ol, fs, "m2.mech.F_z", "servo.p2.o")) < 1.0e-9
end
