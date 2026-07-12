# # Fabry-Perot AC transfer function
#
# The AC analysis route through the ecosystem path (PLAN_rationalize_late §4,
# validation 4): numeric linearization at the operating point
# (`linearize_proper`) → ControlSystemsBase (`ss` / `freqresp` / `poles` /
# `dcgain`). OMTK writes no frequency-sweep code of its own. Two independent
# checks:
#
# 1. the representation-parameterized analytic reference
#    ```math
#    H(s) = \frac{t_1}{1 + r_1 r_2 P(s)^2}, \qquad P(s) = e^{ikL}\,\mathrm{AP}_n(s)
#    ```
#    with ``\mathrm{AP}_n`` the order-``n`` all-pass cascade (here ``n = 1``);
# 2. the hand-rolled symbolic path (`linearize_symbolic` → `fixpoint_sub` →
#    ``C(i\omega I - A)^{-1}B + D`` per frequency), kept as a pedagogical oracle
#    only.
#
# The input is the laser's quadrature modulation pair `(δA, δQ)`; the outputs
# are the quadratures of the right-going intracavity field at M1's front face.

using Test
using OpticalModelingToolkit
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ControlSystemsBase: ss, freqresp, poles, dcgain
using Symbolics
using LinearAlgebra

@testset "Fabry-Perot AC — ecosystem path vs analytic and symbolic oracle" begin
    ## Numerical parameters
    R1 = 0.99
    R2 = 0.999
    r1_v = sqrt(R1)
    r2_v = sqrt(R2)
    t1_v = sqrt(1 - R1)
    λ_v  = 1064e-9
    c_v  = 299792458.0
    ## Macroscopic length ≈ 0.1 m, exactly on resonance with our scattering
    ## convention: exp(2 i k L) = -1 ⇒ L = (2n+1) λ/4.
    n_quarter = round(Int, 0.1 / (λ_v / 4))
    L_v = (2 * n_quarter + 1) * λ_v / 4
    τ_v = L_v / c_v

    ## Model: PadeDelay(1) — one all-pass section per direction, 4 states
    @named laser = LaserModulated(A = 1.0, ϕ = 0.0)
    @named m1    = Mirror(R = R1)
    @named m2    = Mirror(R = R2)
    @named L_fwd = FreeSpace(L = L_v, λ = λ_v, rep = PadeDelay(1))
    @named L_bwd = FreeSpace(L = L_v, λ = λ_v, rep = PadeDelay(1))
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

    inputs  = [laser.δA, laser.δQ]
    outputs = [m1.f_out.a[1], m1.f_out.a[2]]

    ## Ecosystem path: numeric linearization → StateSpace → ControlSystemsBase
    mats, _ = linearize_proper(fp, inputs, outputs)
    sys = ss(mats...)

    ## one all-pass state per quadrature per direction
    @test size(mats.A) == (4, 4)
    ## asymptotically stable around the operating point
    @test all(real(p) < 0 for p in poles(sys))

    ωs = 2π .* (10.0 .^ range(2, 8; length = 25))
    resp = freqresp(sys, ωs)
    ## complex field TF from the in-phase drive: δA is the first input
    Hnum = [resp[1, 1, k] + im * resp[2, 1, k] for k in eachindex(ωs)]

    ## Analytic reference, parameterized by the representation.
    ## P(s) = e^{ikL}(1 - sτ/2)/(1 + sτ/2): all-pass, so |H| is exact at every
    ## frequency and only the phase carries the (ωτ)³/12 section error.
    exp2ikL = cis(2 * 2π * L_v / λ_v)  # ≈ -1 by construction
    Hana = function (ω)
        ap = (1 - im * ω * τ_v / 2) / (1 + im * ω * τ_v / 2)
        t1_v / (1 + r1_v * r2_v * exp2ikL * ap^2)
    end
    for (k, ω) in enumerate(ωs)
        @test Hnum[k] ≈ Hana(ω) rtol = 1e-8 atol = 1e-12
    end

    ## DC gain: the algebraic intracavity buildup t1/(1 - r1 r2), identical in
    ## every representation
    dc = dcgain(sys)
    @test dc[1, 1] ≈ t1_v / (1 - r1_v * r2_v) rtol = 1e-10
    @test Hnum[1] ≈ Hana(ωs[1]) rtol = 1e-10

    ## Symbolic oracle: linearize_symbolic → fixpoint_sub → C(iωI−A)⁻¹B + D.
    ## MTK's symbolic path keeps the torn circulating-field loop variables as
    ## states behind differentiated constraints, so B/D come back with
    ## input-derivative columns [u; u̇] and the realization is only valid for
    ## ω ≠ 0 — the reason `linearize_proper` eliminates the algebraic block
    ## statically instead. As an oracle on the sweep it is exact.
    (; A, B, C, D), ssys = ModelingToolkit.linearize_symbolic(fp, inputs, outputs;
        allow_input_derivatives = true)

    pmap = Dict{Any, Any}()
    for p in parameters(ssys)
        if ModelingToolkit.hasdefault(p)
            pmap[p] = ModelingToolkit.getdefault(p)
        end
    end
    sub = (M) -> Symbolics.value.(ModelingToolkit.fixpoint_sub(M, pmap; fold = Val(true)))
    A_n, B_n, C_n, D_n = Float64.(sub(A)), Float64.(sub(B)), Float64.(sub(C)), Float64.(sub(D))
    @test all(isfinite, A_n) && all(isfinite, B_n) &&
          all(isfinite, C_n) && all(isfinite, D_n)

    nu = length(inputs)
    for (k, ω) in enumerate(ωs)
        v = C_n * (((im * ω) * I - A_n) \ B_n) + D_n
        vu = size(B_n, 2) == 2nu ? v[:, 1:nu] .+ (im * ω) .* v[:, (nu + 1):end] : v
        ## the two paths agree to the finesse-amplified roundoff floor (~1e-8)
        @test Hnum[k] ≈ vu[1, 1] + im * vu[2, 1] rtol = 1e-6 atol = 1e-12
    end

    println("✓ Fabry-Perot AC: ecosystem path matches analytic reference and ",
        "symbolic oracle over ", length(ωs), " points; |H(0)| = ",
        round(abs(dc[1, 1]), digits = 3))
end
