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
using ControlSystemsBase: ss, poles, dcgain
using SciMLBase

const kat_suspended = read(joinpath(pkgdir(OpticalModelingToolkit), "models", "suspended_fp.kat"), String)

@testset "parse_kat" begin
    els = parse_kat(kat_suspended)
    @test [e.kind for e in els] == [:laser, :mirror, :mirror, :space, :space, :pendulum]
    @test els[2].params[:R] == 0.99
    @test els[end].nodes == [(:m2, :mech)]
    @test els[end].params[:Qz] == 2.0
    ## analysis lines are accepted and skipped
    @test length(parse_kat("laser l P=1\nfsig(1)\nfrequency_response(1, a.b, c.d)\n")) == 1
    ## unknown elements and non-numeric values are errors
    @test_throws ArgumentError parse_kat("beamsplitter b R=0.5")
    @test_throws ArgumentError parse_kat("laser l P=1k")
    @test_throws ArgumentError parse_kat("mirror m R=0.9 phi=deg")
end

@testset "build_kat: suspended cavity at resonance" begin
    M, fz, Q = 0.1, 20.0, 2.0
    ω0 = 2π * fz
    plant, h = build_kat(kat_suspended)
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
