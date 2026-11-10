# # KatScript models
#
# `models/*.kat` are the single source of truth. This file checks the subset of
# KatScript that OMT reads: the parser, and that the built network has the
# suspended-cavity dynamics that the file describes. The comparison with Finesse.jl
# on the same file is in `test/finesse/`.

using Test
using OpticalModelingToolkit
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ControlSystemsBase: ss, poles, dcgain, freqresp
using SciMLBase

const kat_suspended = read(joinpath(pkgdir(OpticalModelingToolkit), "models", "suspended_fp.kat"), String)

@testset "parse_kat" begin
    els = parse_kat(kat_suspended)
    @test [e.kind for e in els] == [:laser, :mirror, :mirror, :space, :space, :pendulum]
    @test els[2].params[:R] == 0.99
    @test els[end].nodes == [(:m2, :mech)]
    @test els[end].params[:Qz] == 2.0
    @test els[3].params[:phi] == 0.5
    ## analysis lines are accepted and skipped
    @test length(parse_kat("laser l P=1\nfsig(1)\nfrequency_response(1, a.b, c.d)\n")) == 1
    ## unknown elements and non-numeric values are errors
    @test_throws ArgumentError parse_kat("beamsplitter b R=0.5")
    @test_throws ArgumentError parse_kat("laser l P=1k")
    @test_throws ArgumentError parse_kat("mirror m R=0.9 phi=deg")
    ## lists, with or without spaces, positional or named
    z = only(parse_kat("zpk f [-1, -2.5] [] gain=3"))
    @test z.lists[:_1] == [-1.0, -2.5] && z.lists[:_2] == Float64[] && z.params[:gain] == 3
    @test only(parse_kat("zpk f z=[] p=[-10] k=1")).lists[:p] == [-10.0]
    ## dotted signal nodes keep their full path
    w = only(parse_kat("wire w rd.DC.o m2.mech.F_z"))
    @test w.nodes == [(:rd, Symbol("DC.o")), (:m2, Symbol("mech.F_z"))]
    ## a signal `link` is a chain of wires; an optical one is refused, not dropped
    ws = parse_kat("link(rd.DC.o, f.p1.i)\nlink(f.p2.o, g.p1.i, m2.mech.F_z)")
    @test [e.kind for e in ws] == [:wire, :wire, :wire]
    @test ws[3].nodes == [(:g, Symbol("p1.i")), (:m2, Symbol("mech.F_z"))]
    @test_throws ArgumentError parse_kat("link(l1, m1)")
end

@testset "ZPK matches its transfer function" begin
    s(f) = 2π * im * f
    ## Finesse's default k gives unit DC gain, so `gain` is the DC gain
    for (z, p, gain) in (
            (Float64[], [-10.0], 3.0), ([-188.5], [-6283.2], -2.0),
            ([-1.0, -3.0], [-2.0, -50.0, -400.0], 0.5), ([-5.0], [-7.0], 1.0),
        )
        @named f = ZPK(; z, p, gain)
        @named w = System(Equation[], t; systems = [f])
        mats, _ = linearize_proper(w, [f.u], [f.y])
        H(fr) = gain * prod(abs, p) / prod(abs, z; init = 1.0) *
            prod((s(fr) - zi for zi in z); init = 1.0 + 0im) / prod(s(fr) - pj for pj in p)
        for fr in (0.1, 3.0, 300.0)
            @test only(freqresp(ss(mats...), [2π * fr])) ≈ H(fr) rtol = 1.0e-10
        end
    end
    ## a zero at the origin makes k = 1, as in Finesse
    @named hp = ZPK(; z = [0.0], p = [-10.0])
    @named w = System(Equation[], t; systems = [hp])
    mats, _ = linearize_proper(w, [hp.u], [hp.y])
    @test only(freqresp(ss(mats...), [2π * 1.0])) ≈ s(1.0) / (s(1.0) + 10) rtol = 1.0e-10
    @test_throws ArgumentError ZPK(; name = :bad, z = [-1.0, -2.0], p = [-3.0])
end

@testset "build_kat: suspended cavity at resonance" begin
    M, fz, Q = 0.1, 20.0, 2.0
    ω0 = 2π * fz
    plant, h = build_kat(replace(kat_suspended, "phi=0.5" => "phi=0"))
    ## at zero detuning the optical spring is zero, so the poles are the bare pendulum's
    mats, _ = linearize_proper(plant, [h.F_z], [h.z]; op = Dict(h.z => 0.0, h.v => 0.0))
    expected = let b = ω0 / Q, d = sqrt(complex(b^2 - 4ω0^2))
        [(-b + d) / 2, (-b - d) / 2]
    end
    @test sort(poles(ss(mats...)); by = z -> (imag(z), real(z))) ≈
        sort(expected; by = z -> (imag(z), real(z))) rtol = 1.0e-6
    ## DC: a static force moves a free-spring suspension by F/k
    @test only(dcgain(ss(mats...))) ≈ 1 / (M * ω0^2) rtol = 1.0e-6
end
