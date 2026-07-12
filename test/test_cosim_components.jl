# # Co-sim components: Photodetector and MovableMirror
#
# The same on-resonance Fabry-Perot as [the DC operating point](@ref
# "Fabry-Perot DC operating point"), but with the end mirror movable and a
# photodetector on transmission. At `x = 0` the DC solution must reproduce the
# analytic Fabry-Perot transmission exactly; sweeping `x` must trace the Airy
# fringe
#
# ```math
# P_t(x) = \left| \frac{t_1 t_2 A}{1 + r_1 r_2 e^{2ikL} e^{2ikx}} \right|^2 .
# ```

using Test
using OpticalModelingToolkit
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using NonlinearSolve
using SciMLBase

@testset "co-sim components — movable-mirror Fabry-Perot with PD" begin
    R1 = 0.99
    R2 = 0.999
    r1 = sqrt(R1)
    r2 = sqrt(R2)
    t1 = sqrt(1 - R1)
    t2 = sqrt(1 - R2)
    λ = 1064.0e-9
    A = 1.0
    ## With the S = [[t, r], [-r, t]] convention the loop gain is
    ## -r1 r2 e^{2ik(L+x)}, so L = λ/2 (e^{2ikL} = +1) sits at the transmission
    ## *minimum*; the test checks formula equality, not resonance.
    L = λ / 2

    @named laser = Laser(A = A, ϕ = 0.0)
    @named m1 = Mirror(R = R1)
    @named m2 = MovableMirror(R = R2, λ = λ)
    @named L_fwd = FreeSpace(L = L, λ = λ)
    @named L_bwd = FreeSpace(L = L, λ = λ)
    @named pd = Photodetector()

    ## analytic transmitted power vs microscopic end-mirror position x:
    ## circulating (right-going at M2 input) E = t1 A e^{ikL} / (1 - r1 r2 e^{2ikL} e^{2ikx} · (-1)²)
    ## with the S = [[t, r], [-r, t]] convention both reflections inside the
    ## cavity contribute one (-r1)·(r2) → net minus sign:
    ## E_circ = t1 A e^{ikL} / (1 + r1 r2 e^{2i(kL + kx)}), P_t = |t2 E_circ|²
    P_analytic(x) = begin
        φ = 2 * (2π / λ) * (L + x)
        E = t1 * A / abs(1 + r1 * r2 * cis(φ))
        (t2 * E)^2
    end

    function transmitted_power(xval)
        eqs = [
            connect(laser.port_out, m1.b_in),
            connect(m1.b_out, laser.port_in),
            connect(m1.f_out, L_fwd.port_in),
            connect(L_fwd.port_out, m2.b_in),
            connect(m2.b_out, L_bwd.port_in),
            connect(L_bwd.port_out, m1.f_in),
            connect(m2.f_out, pd.port_in),
            connect(pd.port_out, m2.f_in),
            m2.x ~ xval,
        ]
        @named fp = System(eqs, t; systems = [laser, m1, m2, L_fwd, L_bwd, pd])
        sys = mtkcompile(fp)
        if isempty(unknowns(sys))
            ## fully torn: evaluate the observed PD power directly
            prob = NonlinearProblem(sys, Dict())
            sol = solve(prob)
            return sol[pd.P]
        else
            prob = NonlinearProblem(sys, Dict(u => 0.0 for u in unknowns(sys)))
            sol = solve(prob)
            @test SciMLBase.successful_retcode(sol)
            return sol[pd.P]
        end
    end

    ## at x = 0: MovableMirror must reduce exactly to Mirror
    @test transmitted_power(0.0) ≈ P_analytic(0.0) rtol = 1.0e-9

    ## a handful of microscopic offsets
    for x in (0.05e-9, -0.1e-9, 0.25e-9, 1.0e-9)
        @test transmitted_power(x) ≈ P_analytic(x) rtol = 1.0e-6
    end
end
