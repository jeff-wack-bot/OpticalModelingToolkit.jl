# Phase 1a — Fabry-Perot DC operating point.
#
# Wires up the simplest meaningful network exercising directed optical
# connectors: two mirrors face-to-face with explicit forward and backward
# Length elements for free-space propagation, a Laser on the input side, and
# a Vacuum termination on the output side.
#
#     Laser ──> M1 ──L_fwd──> M2 ──> Vacuum
#                ^            |
#                └───L_bwd────┘
#
# The whole system is purely algebraic, so `mtkcompile` produces an
# equation-only system that we solve as a `NonlinearProblem`. The intracavity
# field amplitude must match the analytical Fabry-Perot formula
#
#     E_cav = t1 * A / (1 + r1 * r2 * exp(2 i k L))
#
# at machine precision.

using Test
using MTKOptical
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using NonlinearSolve
using SciMLBase

@testset "Phase 1a — Fabry-Perot DC operating point" begin
    # ----------------------------------------------------------------------
    # Numerical parameters
    # ----------------------------------------------------------------------
    R1 = 0.99
    R2 = 0.999
    r1 = sqrt(R1)
    r2 = sqrt(R2)
    t1 = sqrt(1 - R1)
    t2 = sqrt(1 - R2)
    λ = 1064e-9
    A = 1.0

    # Pick a length that puts the cavity exactly on resonance, so that
    # exp(2 i k L) = -1 ⇒ E_cav = t1 * A / (1 - r1 * r2). Any half-integer
    # multiple of λ/2 works; we use a single half-wavelength because the
    # algebraic Length component is independent of L modulo λ/2 anyway.
    L = λ / 2

    # ----------------------------------------------------------------------
    # Build the model
    # ----------------------------------------------------------------------
    @named laser = Laser(A = A, ϕ = 0.0)
    @named m1    = Mirror(R = R1)
    @named m2    = Mirror(R = R2)
    @named L_fwd = Length(L = L, λ = λ)
    @named L_bwd = Length(L = L, λ = λ)
    @named vac   = Vacuum()

    eqs = [
        # Laser into back of M1, M1 reflection into laser sink
        connect(laser.port_out, m1.b_in),
        connect(m1.b_out,       laser.port_in),

        # Forward propagation: M1 front out → M2 back in
        connect(m1.f_out,       L_fwd.port_in),
        connect(L_fwd.port_out, m2.b_in),

        # Backward propagation: M2 back out → M1 front in
        connect(m2.b_out,       L_bwd.port_in),
        connect(L_bwd.port_out, m1.f_in),

        # Vacuum on the far side of M2
        connect(m2.f_out, vac.port_in),
        connect(vac.port_out, m2.f_in),
    ]

    @named fp = System(eqs, t; systems = [laser, m1, m2, L_fwd, L_bwd, vac])
    sys = mtkcompile(fp)

    # ----------------------------------------------------------------------
    # Solve as a NonlinearProblem
    # ----------------------------------------------------------------------
    # The directed-connector graph leaves a handful of array-component
    # unknowns after `mtkcompile`. Zero is a perfectly good initial guess for
    # an algebraic linear system, so just hand back zeros for whatever
    # `unknowns(sys)` reports as still free.
    u0 = [u => 0.0 for u in unknowns(sys)]
    prob = NonlinearProblem(sys, u0)
    sol = solve(prob, NewtonRaphson(); reltol = 1e-14, abstol = 1e-14)
    @test SciMLBase.successful_retcode(sol)

    # ----------------------------------------------------------------------
    # Analytical reference
    # ----------------------------------------------------------------------
    # Right-going amplitude just inside the cavity (i.e. at m1.f_out).
    E1plus_analytic = t1 * A / (1 + r1 * r2 * cis(2 * 2π * L / λ))

    # Reflected amplitude back toward the laser (m1.b_out).
    E1minus_analytic = r2 * cis(2 * 2π * L / λ) * E1plus_analytic
    E_refl_analytic  = t1 * E1minus_analytic + r1 * A

    # Transmitted amplitude (m2.f_out).
    E_trans_analytic = t2 * cis(2π * L / λ) * E1plus_analytic

    # ----------------------------------------------------------------------
    # Compare
    # ----------------------------------------------------------------------
    E1plus_re = sol[m1.f_out.a[1]]
    E1plus_im = sol[m1.f_out.a[2]]
    @test E1plus_re ≈ real(E1plus_analytic) rtol=1e-10
    @test E1plus_im ≈ imag(E1plus_analytic) rtol=1e-10 atol=1e-12

    E_refl_re = sol[m1.b_out.a[1]]
    E_refl_im = sol[m1.b_out.a[2]]
    @test E_refl_re ≈ real(E_refl_analytic) rtol=1e-10
    @test E_refl_im ≈ imag(E_refl_analytic) rtol=1e-10 atol=1e-12

    E_trans_re = sol[m2.f_out.a[1]]
    E_trans_im = sol[m2.f_out.a[2]]
    @test E_trans_re ≈ real(E_trans_analytic) rtol=1e-10
    @test E_trans_im ≈ imag(E_trans_analytic) rtol=1e-10 atol=1e-12

    # Energy conservation: |E_refl|² + |E_trans|² ≈ |A|²
    P_in    = A^2
    P_refl  = E_refl_re^2 + E_refl_im^2
    P_trans = E_trans_re^2 + E_trans_im^2
    @test P_refl + P_trans ≈ P_in rtol=1e-10

    println("✓ Fabry-Perot DC: |E_cav|² = ", round(E1plus_re^2 + E1plus_im^2, digits=4),
            "  P_refl = ", round(P_refl, digits=6),
            "  P_trans = ", round(P_trans, digits=6))
end
