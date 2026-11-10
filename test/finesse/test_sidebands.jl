# # Exact sideband transfer functions
#
# `to_scattering` keeps every `FreeSpace` as an exact delay, so its `G(2πif)`
# is the same object as Finesse's upper-sideband optical transfer function
# (`freqresp3`) at `fsig.f = f`. Both are transcendental in `f`, so the
# comparison runs through the cavity pole, out past the free spectral range,
# and onto FSR multiples, where a wrong delay or a wrong FSR shows at once.
#
# The input is the field entering the first optic from the laser. OMT names it
# after the source terminal; Finesse names it by the node it enters.
#
# ## The one modeling difference: tuning dispersion
#
# Finesse treats a tuning as a displacement, so the sideband at `f₀ + f`
# reflects with phase `2φ(1 + f/f₀)`. OMT applies the carrier's `2ϕ` to every
# sideband. For a 1° tuning at 1 GHz that is a 6e-8 rad difference per
# reflection, multiplied by the cavity's finesse. The detuned comparisons
# account for it exactly: at each frequency the OMT model is rebuilt with
# every tuning scaled by `1 + f/f₀`, which makes the two codes compute the same
# thing. Undetuned networks need no correction (`0 · (1 + f/f₀) = 0`), and
# the last testset shows the correction is not cosmetic.

const C0 = 299792458.0

disp(ϕ, f) = ϕ * (1 + f / F0)

@testset "Fabry-Perot reflection and transmission" begin
    for (ϕ1, ϕ2, L) in ((0.0, 0.0, 0.1), (0.0, deg2rad(1.0), 0.1), (0.3, 0.4, 3.0))
        p = (; R1 = 0.99, T1 = 0.0099, R2 = 0.999, T2 = 0.0008, ϕ1, ϕ2, P = 1.0, L)
        fsr = C0 / 2L
        for f in (1.0, 1.0e3, 1.0e5, 1.0e6, 3.7e6, 0.3fsr, 0.5fsr, fsr, fsr + 1.0e5, 2.25fsr)
            fp, _ = omt_fp(; p..., ϕ1 = disp(ϕ1, f), ϕ2 = disp(ϕ2, f))
            parts = to_scattering(fp)
            iu = findfirst(==(:l1), parts.u)
            irefl = findfirst(==(:l1), parts.y)
            itrans = findfirst(==(:pd), parts.y)
            G = omt_tf(parts, f)
            fin = finesse_freqresp3(
                kat_fp(; p...), f, "$(mport("m1", 2)).i",
                ["$(mport("m1", 2)).o", "$(mport("m2", 2)).o"]
            )
            @test G[irefl, iu] ≈ fin[1] rtol = 1.0e-9
            @test G[itrans, iu] ≈ fin[2] rtol = 1.0e-9
        end
    end
end

@testset "Michelson antisymmetric and symmetric ports" begin
    kw = (; R = 0.4, ϕx = 0.2, ϕy = -0.7, ϕbs = 0.1, Rx = 0.99, Ry = 0.95)
    src = """
    $(kat_bs("bs"; R = kw.R, ϕ = kw.ϕbs))
    $(kat_mirror("ex"; R = kw.Rx, ϕ = kw.ϕx))
    $(kat_mirror("ey"; R = kw.Ry, ϕ = kw.ϕy))
    space sy $(bsport("bs", 2)) $(mport("ey", 1)) L=1.3
    space sx $(bsport("bs", 3)) $(mport("ex", 1)) L=1
    """
    for f in (1.0e3, 1.0e6, 7.7e7, 1.0e8, 3.0e8)
        mich, _ = omt_michelson(;
            kw..., ϕx = disp(kw.ϕx, f), ϕy = disp(kw.ϕy, f), ϕbs = disp(kw.ϕbs, f)
        )
        parts = to_scattering(mich)
        iu = findfirst(==(:l1), parts.u)
        ias = findfirst(==(:pd), parts.y)
        isym = findfirst(==(:l1), parts.y)
        G = omt_tf(parts, f)
        fin = finesse_freqresp3(
            src, f, "$(bsport("bs", 1)).i", ["$(bsport("bs", 4)).o", "$(bsport("bs", 1)).o"]
        )
        @test G[ias, iu] ≈ fin[1] rtol = 1.0e-10 atol = 1.0e-14
        @test G[isym, iu] ≈ fin[2] rtol = 1.0e-10 atol = 1.0e-14
    end
end

@testset "coupled-cavity reflection and transmission" begin
    for (ϕ2, ϕ3) in ((0.0, 0.0), (0.02, -0.01))
        for f in (1.0e2, 1.0e5, 1.0e6, 3.0e7, 1.5e8, 4.9e8)
            cc, _ = omt_coupled(; ϕ2 = disp(ϕ2, f), ϕ3 = disp(ϕ3, f))
            parts = to_scattering(cc)
            iu = findfirst(==(:l1), parts.u)
            it = findfirst(==(:pd), parts.y)
            ir = findfirst(==(:l1), parts.y)
            G = omt_tf(parts, f)
            fin = finesse_freqresp3(
                kat_coupled(; ϕ2, ϕ3), f, "$(mport("m1", 2)).i",
                ["$(mport("m3", 2)).o", "$(mport("m1", 2)).o"]
            )
            @test G[it, iu] ≈ fin[1] rtol = 1.0e-9
            @test G[ir, iu] ≈ fin[2] rtol = 1.0e-9
        end
    end
end

# Without the correction, a detuned high-finesse cavity near the FSR differs
# by far more than the tolerances above: the dispersion is real physics that
# OMT's tuning model leaves out, and these tests would catch either code
# changing its mind about it.

@testset "tuning dispersion is the whole difference" begin
    p = (; R1 = 0.99, T1 = 0.0099, R2 = 0.999, T2 = 0.0008, ϕ1 = 0.0, ϕ2 = deg2rad(1.0), P = 1.0, L = 0.1)
    f = C0 / 2p.L + 1.0e5
    fin = finesse_freqresp3(kat_fp(; p...), f, "$(mport("m1", 2)).i", ["$(mport("m2", 2)).o"])[1]
    fp, _ = omt_fp(; p...)
    parts = to_scattering(fp)
    G = omt_tf(parts, f)[findfirst(==(:pd), parts.y), findfirst(==(:l1), parts.u)]
    @test !isapprox(G, fin; rtol = 1.0e-6)
    @test isapprox(G, fin; rtol = 1.0e-3)
end
