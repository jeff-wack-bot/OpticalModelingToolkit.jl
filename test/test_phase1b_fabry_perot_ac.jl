# Phase 1b — Fabry-Perot AC transfer function via linearize_symbolic.
#
# Same Fabry-Perot topology as Phase 1a, but the algebraic `Length` elements
# are replaced by `LengthDelay`, which models free-space propagation as a
# first-order single-pole approximation
#
#     τ ȯ + o = exp(ikL) i,    τ = L / c
#
# The cavity now has 4 state variables (2 quadratures × 2 length elements),
# and `linearize_symbolic` produces a non-trivial `(A, B, C, D)`. The input
# is the small-signal amplitude perturbation `δA(t)` of `LaserModulated`; the
# outputs are the in-phase and quadrature components of the right-going field
# at M1's front face (i.e. the intracavity right-going amplitude).
#
# The transfer function from `δA` to `m1.f_out.a` is, analytically,
#
#     H(ω) = t1 (1 + iωτ)² / [(1 + iωτ)² + r1 r2 exp(2 i k L)]
#
# evaluated at the carrier phase. We compare this against
#
#     H_num(ω) = C (iωI − A)⁻¹ B + D
#
# at every point of a logarithmic frequency sweep.

using Test
using OpticalModelingToolkit
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using NonlinearSolve
using SciMLBase
using Symbolics
using LinearAlgebra

@testset "Phase 1b — Fabry-Perot AC transfer function" begin
    # ----------------------------------------------------------------------
    # Numerical parameters
    # ----------------------------------------------------------------------
    R1 = 0.99
    R2 = 0.999
    r1_v = sqrt(R1)
    r2_v = sqrt(R2)
    t1_v = sqrt(1 - R1)
    t2_v = sqrt(1 - R2)
    λ_v  = 1064e-9
    A_v  = 1.0
    c_v  = 299792458.0
    # Macroscopic cavity length so that c/L is well above the audio band
    # (cavity pole around c/(L) ≈ 3 GHz). Pick a length that's also exactly
    # on resonance with our scattering convention: exp(2 i k L) = -1 ⇒ kL =
    # π/2 + nπ ⇒ L = (2n+1) λ/4. We use n large enough to be physically
    # plausible.
    n_quarter = round(Int, 0.1 / (λ_v / 4))     # ≈ 0.1 m
    L_v = (2 * n_quarter + 1) * λ_v / 4

    # ----------------------------------------------------------------------
    # Build the model with the dynamic LengthDelay
    # ----------------------------------------------------------------------
    @named laser = LaserModulated(A = A_v, ϕ = 0.0)
    @named m1    = Mirror(R = R1)
    @named m2    = Mirror(R = R2)
    @named L_fwd = LengthDelay(L = L_v, λ = λ_v, c = c_v)
    @named L_bwd = LengthDelay(L = L_v, λ = λ_v, c = c_v)
    @named vac   = Vacuum()

    eqs = [
        connect(laser.port_out, m1.b_in),
        connect(m1.b_out,       laser.port_in),

        connect(m1.f_out,       L_fwd.port_in),
        connect(L_fwd.port_out, m2.b_in),

        connect(m2.b_out,       L_bwd.port_in),
        connect(L_bwd.port_out, m1.f_in),

        connect(m2.f_out, vac.port_in),
        connect(vac.port_out, m2.f_in),
    ]

    @named fp = System(eqs, t; systems = [laser, m1, m2, L_fwd, L_bwd, vac])

    # ----------------------------------------------------------------------
    # Symbolic linearization
    # ----------------------------------------------------------------------
    # Inputs: laser amplitude perturbation. Outputs: intracavity right-going
    # field at M1's front face, in-phase and quadrature.
    inputs  = [laser.δA]
    outputs = [m1.f_out.a[1], m1.f_out.a[2]]

    (; A, B, C, D), ssys = ModelingToolkit.linearize_symbolic(fp, inputs, outputs)

    # The system has two state variables per LengthDelay quadrature, so we
    # expect 4 states, 1 input, 2 outputs.
    @test size(A) == (4, 4)
    @test size(B) == (4, 1)
    @test size(C) == (2, 4)
    @test size(D) == (2, 1)

    # ----------------------------------------------------------------------
    # Substitute parameter values to get numeric matrices
    # ----------------------------------------------------------------------
    # Build a substitution map from every parameter in the simplified system
    # to its constructor-time default. The Mirror, LengthDelay, and
    # LaserModulated components all bake their numeric R/L/λ/c/A values into
    # parameter defaults at construction time, so this is a one-shot
    # substitution.
    pmap = Dict{Any, Any}()
    for p in parameters(ssys)
        if ModelingToolkit.hasdefault(p)
            pmap[p] = ModelingToolkit.getdefault(p)
        end
    end
    sub = (M) -> Symbolics.value.(ModelingToolkit.fixpoint_sub(M, pmap; fold = Val(true)))

    A_n = Float64.(sub(A))
    B_n = Float64.(sub(B))
    C_n = Float64.(sub(C))
    D_n = Float64.(sub(D))

    # All four matrices must come out as concrete numbers — anything left
    # symbolic indicates a parameter we forgot to bind.
    @test all(isfinite, A_n)
    @test all(isfinite, B_n)
    @test all(isfinite, C_n)
    @test all(isfinite, D_n)

    # The system should be asymptotically stable around the operating point:
    # all eigenvalues of A in the open left half-plane.
    eigs = eigvals(A_n)
    @test all(real(e) < 0 for e in eigs)

    # ----------------------------------------------------------------------
    # Numerical transfer function H(iω) = C (iωI − A)⁻¹ B + D
    # ----------------------------------------------------------------------
    Hnum = function (ω)
        Y = (im * ω) * I - A_n
        v = C_n * (Y \ B_n) + D_n
        # Output is a 2-vector whose entries are themselves complex
        # frequency responses: v[1] = response of the real quadrature, v[2]
        # = response of the imaginary quadrature. Combine them into a single
        # complex field transfer function.
        v[1, 1] + im * v[2, 1]
    end

    # ----------------------------------------------------------------------
    # Analytical reference
    # ----------------------------------------------------------------------
    # The single-pole length approximation has the same s → 0 limit as the
    # algebraic case, so the analytical transfer function is
    #
    #   H(s) = t1 (1 + sτ)² / [(1 + sτ)² + r1 r2 exp(2ikL)]
    #
    # evaluated at s = iω. With our resonance tuning, exp(2ikL) = -1.
    τ_v = L_v / c_v
    exp2ikL = cis(2 * 2π * L_v / λ_v)  # ≈ -1 by construction
    Hana = function (ω)
        x = 1 + im * ω * τ_v
        t1_v * x^2 / (x^2 + r1_v * r2_v * exp2ikL)
    end

    # ----------------------------------------------------------------------
    # Sweep
    # ----------------------------------------------------------------------
    # Sweep from 100 Hz to 100 MHz (well below the cavity FSR ≈ c/(2L) ≈
    # 1.5 GHz, so the single-pole approximation is reasonable). Compare
    # numeric vs. analytic at every point.
    ωs = 2π .* (10.0 .^ range(2, 8; length = 25))
    max_relerr = 0.0
    for ω in ωs
        h_num = Hnum(ω)
        h_ana = Hana(ω)
        relerr = abs(h_num - h_ana) / abs(h_ana)
        max_relerr = max(max_relerr, relerr)
        @test isapprox(h_num, h_ana; rtol = 1e-8, atol = 1e-12)
    end

    # DC value sanity check: should equal the algebraic intracavity buildup
    # (linear gain from δA to the field at DC is t1 / (1 - r1 r2)). The
    # imaginary part is not exactly zero because `2π L / λ` in float
    # arithmetic isn't exactly an odd multiple of π — the analytical formula
    # has the same residual, so we just compare against `Hana(0)`.
    h_dc_num = Hnum(0.0)
    h_dc_ana = Hana(0.0)
    @test real(h_dc_num) ≈ t1_v / (1 - r1_v * r2_v) rtol=1e-10
    @test h_dc_num ≈ h_dc_ana rtol=1e-10 atol=1e-12

    println("✓ Fabry-Perot AC: states=", size(A_n, 1),
            "  |H(0)|=", round(abs(h_dc_num), digits=3),
            "  max relerr over sweep=", round(max_relerr, sigdigits=3))
end
