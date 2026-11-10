# # DC fields in composed networks
#
# The same networks built in both codes, solved at DC, and compared field by
# field. OMT solves the realified MTK equations (`mtkcompile` +
# `NonlinearProblem`); Finesse solves its sparse carrier matrix. Nothing is
# analytic here: the reference is the other simulator.

# ## Fabry-Perot, lossy and detuned
#
# ```
#     l1 ── p2 [m1] p1 ══ L ══ p1 [m2] p2 ── pd
# ```

function omt_fp(; R1, T1, R2, T2, ϕ1, ϕ2, P, L)
    @named l1 = Laser(A = sqrt(P))
    @named m1 = Mirror(R = R1, T = T1, ϕ = ϕ1)
    @named sp = FreeSpace(; L)
    @named m2 = Mirror(R = R2, T = T2, ϕ = ϕ2)
    @named pd = Photodetector()
    eqs = vcat(link(l1.p1, m1.p2), link(m1.p1, sp.p1), link(sp.p2, m2.p1), link(m2.p2, pd.p1))
    return System(eqs, t; systems = [l1, m1, sp, m2, pd], name = :fp), (; m1, m2, pd)
end

kat_fp(; R1, T1, R2, T2, ϕ1, ϕ2, P, L) = """
laser l1 P=$P
$(kat_mirror("m1"; R = R1, T = T1, ϕ = ϕ1))
$(kat_mirror("m2"; R = R2, T = T2, ϕ = ϕ2))
space s0 l1.p1 $(mport("m1", 2)) L=1
space sc $(mport("m1", 1)) $(mport("m2", 1)) L=$L
"""

@testset "Fabry-Perot fields at every port" begin
    base = (; R1 = 0.99, T1 = 0.0099, R2 = 0.999, T2 = 0.0008, ϕ1 = 0.0, P = 2.0, L = 0.1)
    for ϕ2 in (0.0, 1.0e-3, -4.0e-3, 0.05, π / 4, π / 2), ϕ1 in (0.0, 0.3)
        p = (; base..., ϕ1, ϕ2)
        fp, c = omt_fp(; p...)
        sol = omt_dc(fp)
        fin = finesse_run(
            kat_fp(; p...) * """
                ad circ $(mport("m1", 1)).o 0
                ad refl $(mport("m1", 2)).o 0
                ad atm2 $(mport("m2", 1)).i 0
                ad trans $(mport("m2", 2)).o 0
                pd P $(mport("m2", 2)).o
                noxaxis()
                """
        )
        @test cfield(sol, c.m1.p1.o) ≈ reading(fin, "circ") rtol = 1.0e-10
        @test cfield(sol, c.m1.p2.o) ≈ reading(fin, "refl") rtol = 1.0e-10
        @test cfield(sol, c.m2.p1.i) ≈ reading(fin, "atm2") rtol = 1.0e-10
        @test cfield(sol, c.m2.p2.o) ≈ reading(fin, "trans") rtol = 1.0e-10
        @test sol[c.pd.P] ≈ real(reading(fin, "P")) rtol = 1.0e-10
    end
end

# ## Michelson
#
# A beamsplitter of reflectivity `R` with tuned end mirrors; the detector on
# `p4` is the antisymmetric port, and `p1` (back toward the laser) the
# symmetric one.

function omt_michelson(; R, ϕx, ϕy, ϕbs, Rx, Ry)
    @named l1 = Laser()
    @named bs = Beamsplitter(; R, ϕ = ϕbs)
    @named sx = FreeSpace(L = 1.0)
    @named sy = FreeSpace(L = 1.3)
    @named ex = Mirror(R = Rx, ϕ = ϕx)
    @named ey = Mirror(R = Ry, ϕ = ϕy)
    @named vx = Vacuum()
    @named vy = Vacuum()
    @named pd = Photodetector()
    eqs = vcat(
        link(l1.p1, bs.p1), link(bs.p2, sy.p1), link(sy.p2, ey.p1), link(ey.p2, vy.p1),
        link(bs.p3, sx.p1), link(sx.p2, ex.p1), link(ex.p2, vx.p1), link(bs.p4, pd.p1),
    )
    sys = System(eqs, t; systems = [l1, bs, sx, sy, ex, ey, vx, vy, pd], name = :mich)
    return sys, (; bs, ex, ey, pd)
end

@testset "Michelson fields at both output ports and in the arms" begin
    for (R, ϕx, ϕy, ϕbs, Rx, Ry) in (
            (0.5, 0.0, 0.0, 0.0, 1.0, 1.0), (0.5, 0.0, π / 4, 0.0, 1.0, 1.0),
            (0.4, 0.2, -0.7, 0.1, 0.99, 0.95), (0.3, 1.0, 1.2, -0.4, 0.9, 1.0),
        )
        mich, c = omt_michelson(; R, ϕx, ϕy, ϕbs, Rx, Ry)
        sol = omt_dc(mich)
        fin = finesse_run(
            """
            laser l1 P=1
            $(kat_bs("bs"; R, ϕ = ϕbs))
            $(kat_mirror("ex"; R = Rx, ϕ = ϕx))
            $(kat_mirror("ey"; R = Ry, ϕ = ϕy))
            space s0 l1.p1 $(bsport("bs", 1)) L=1
            space sy $(bsport("bs", 2)) $(mport("ey", 1)) L=1.3
            space sx $(bsport("bs", 3)) $(mport("ex", 1)) L=1
            ad sym $(bsport("bs", 1)).o 0
            ad as $(bsport("bs", 4)).o 0
            ad ax $(mport("ex", 1)).i 0
            ad ay $(mport("ey", 1)).i 0
            pd P $(bsport("bs", 4)).o
            noxaxis()
            """
        )
        @test cfield(sol, c.bs.p1.o) ≈ reading(fin, "sym") rtol = 1.0e-12 atol = 1.0e-14
        @test cfield(sol, c.bs.p4.o) ≈ reading(fin, "as") rtol = 1.0e-12 atol = 1.0e-14
        @test cfield(sol, c.ex.p1.i) ≈ reading(fin, "ax") rtol = 1.0e-12 atol = 1.0e-14
        @test cfield(sol, c.ey.p1.i) ≈ reading(fin, "ay") rtol = 1.0e-12 atol = 1.0e-14
        @test sol[c.pd.P] ≈ real(reading(fin, "P")) rtol = 1.0e-12 atol = 1.0e-14
    end
end

# ## Coupled cavities
#
# Three mirrors in a line: two cavities sharing `m2`. This is the smallest
# network where a field crosses two resonators, and where getting a sign of the
# middle mirror wrong (it is entered from both faces) changes the answer.
#
# ```
#     l1 ── p2 [m1] p1 ══ L1 ══ p1 [m2] p2 ══ L2 ══ p1 [m3] p2 ── pd
# ```

function omt_coupled(; ϕ2, ϕ3)
    @named l1 = Laser()
    @named m1 = Mirror(R = 0.9)
    @named s1 = FreeSpace(L = 1.0)
    @named m2 = Mirror(R = 0.95, T = 0.049, ϕ = ϕ2)
    @named s2 = FreeSpace(L = 0.3)
    @named m3 = Mirror(R = 0.99, ϕ = ϕ3)
    @named pd = Photodetector()
    eqs = vcat(
        link(l1.p1, m1.p2), link(m1.p1, s1.p1), link(s1.p2, m2.p1),
        link(m2.p2, s2.p1), link(s2.p2, m3.p1), link(m3.p2, pd.p1),
    )
    sys = System(eqs, t; systems = [l1, m1, s1, m2, s2, m3, pd], name = :cc)
    return sys, (; m1, m2, m3, pd)
end

kat_coupled(; ϕ2, ϕ3) = """
laser l1 P=1
$(kat_mirror("m1"; R = 0.9))
$(kat_mirror("m2"; R = 0.95, T = 0.049, ϕ = ϕ2))
$(kat_mirror("m3"; R = 0.99, ϕ = ϕ3))
space s0 l1.p1 $(mport("m1", 2)) L=1
space s1 $(mport("m1", 1)) $(mport("m2", 1)) L=1
space s2 $(mport("m2", 2)) $(mport("m3", 1)) L=0.3
"""

@testset "coupled-cavity fields" begin
    for (ϕ2, ϕ3) in ((0.0, 0.0), (0.0, π / 2), (0.02, -0.01), (π / 2, 0.3))
        cc, c = omt_coupled(; ϕ2, ϕ3)
        sol = omt_dc(cc)
        fin = finesse_run(
            kat_coupled(; ϕ2, ϕ3) * """
                ad refl $(mport("m1", 2)).o 0
                ad c1 $(mport("m1", 1)).o 0
                ad c2 $(mport("m2", 2)).o 0
                ad trans $(mport("m3", 2)).o 0
                noxaxis()
                """
        )
        @test cfield(sol, c.m1.p2.o) ≈ reading(fin, "refl") rtol = 1.0e-10
        @test cfield(sol, c.m1.p1.o) ≈ reading(fin, "c1") rtol = 1.0e-10
        @test cfield(sol, c.m2.p2.o) ≈ reading(fin, "c2") rtol = 1.0e-10
        @test cfield(sol, c.m3.p2.o) ≈ reading(fin, "trans") rtol = 1.0e-10
    end
end
