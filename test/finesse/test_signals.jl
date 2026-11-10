# # Linearized signals: the time-domain model's small-signal response
#
# The sideband tests check `to_scattering`, which never compiles the MTK
# model. These check the other path, the one the package exists for: the
# realified DAE, linearized by `linearize_proper`. Finesse computes the same
# responses from its signal matrix.

# A detuned Fabry-Perot with a modulated laser and a photodetector. `m2` is
# movable so that its position can be an input; with `x` bound it is a mirror.

function omt_fp_signal(; rep = Algebraic(), ϕ2, P = 2.0, x = 0.0, laser = Laser)
    @named l1 = laser(A = sqrt(P))
    @named m1 = Mirror(R = 0.99)
    @named sp = FreeSpace(L = 0.1; rep)
    @named m2 = MovableMirror(R = 0.999, ϕ = ϕ2)
    @named pd = Photodetector()
    eqs = vcat(link(l1.p1, m1.p2), link(m1.p1, sp.p1), link(sp.p2, m2.p1), link(m2.p2, pd.p1))
    x === nothing || push!(eqs, m2.x ~ x)
    return System(eqs, t; systems = [l1, m1, sp, m2, pd], name = :fp), (; l1, m2, pd)
end

kat_fp_signal(; ϕ2, P = 2.0) = """
laser l1 P=$P
$(kat_mirror("m1"; R = 0.99))
$(kat_mirror("m2"; R = 0.999, ϕ = ϕ2))
space s0 l1.p1 $(mport("m1", 2)) L=0
space sc $(mport("m1", 1)) $(mport("m2", 1)) L=0.1
"""

# ## Laser amplitude modulation to detected power
#
# Finesse's `amp` drive is relative amplitude modulation `δA/A`, and `pd1` at
# `fsig` reads the amplitude of the resulting power oscillation. OMT's
# `LaserModulated.δA` is absolute, so the Finesse response is OMT's
# `∂P/∂δA` times `A`. Both quadrature inputs are named, because an unbound
# `[input = true]` variable makes `mtkcompile` fail.

finesse_am(ϕ2, f) = reading(
    finesse_run(
        kat_fp_signal(; ϕ2) * """
            fsig($f)
            sgen sg l1.amp 1 0
            pd1 P $(mport("m2", 2)).o fsig
            noxaxis()
            """
    ), "P"
)

@testset "laser AM to transmitted power" begin
    A = sqrt(2.0)
    for ϕ2 in (0.0, deg2rad(1.0), deg2rad(-3.0))
        ## Algebraic propagation: the DC gain, exactly
        fp, c = omt_fp_signal(; ϕ2, laser = LaserModulated)
        mats, _ = linearize_proper(fp, [c.l1.δA, c.l1.δQ], [c.pd.P])
        @test dcgain(ss(mats...))[1, 1] * A ≈ finesse_am(ϕ2, 1.0e-3) rtol = 1.0e-9

        ## PadeDelay(8): the full frequency response through the cavity pole
        ## (~1.2 MHz) and past it; the all-pass phase error at 10 MHz is
        ## (2πfτ)³/(12·8²) ~ 1e-8 per pass, and the cavity's finesse
        ## multiplies that by ~600.
        fp, c = omt_fp_signal(; ϕ2, laser = LaserModulated, rep = PadeDelay(8))
        mats, _ = linearize_proper(fp, [c.l1.δA, c.l1.δQ], [c.pd.P])
        sys = ss(mats...)
        for f in (1.0e3, 1.0e5, 1.0e6, 3.0e6, 1.0e7)
            @test freqresp(sys, 2π * f)[1, 1, 1] * A ≈ finesse_am(ϕ2, f) rtol = 1.0e-5
        end
    end
end

# ## The optical spring
#
# Moving the detuned end mirror changes the circulating power and so the
# radiation-pressure force. OMT linearizes `m2.x → m2.F`; Finesse reports the
# response from `m2.mech.z` to `m2.mech.F_z`. The two agree in magnitude and
# have opposite sign: with both codes' sign conventions as they stand, Finesse's
# `z → F_z` is `-∂F/∂x` of OMT's `MovableMirror`, i.e. the stiffness `K`.

@testset "optical spring constant" begin
    for ϕ2 in (deg2rad(1.0), deg2rad(-0.2), deg2rad(0.05), 0.0, deg2rad(5.0))
        fp, c = omt_fp_signal(; ϕ2, x = nothing)
        mats, _ = linearize_proper(fp, [c.m2.x], [c.m2.F])
        dFdx = dcgain(ss(mats...))[1, 1]
        sol = finesse_run(kat_fp_signal(; ϕ2) * "\nfsig(1)\nfreqresp(1e-3, [m2.mech.z], [m2.mech.F_z])\n")
        K = only(sol.data)
        @test -dFdx ≈ real(K) rtol = 1.0e-8 atol = 1.0e-6
        @test abs(imag(K)) < 1.0e-6 * max(abs(real(K)), 1)
    end
end
