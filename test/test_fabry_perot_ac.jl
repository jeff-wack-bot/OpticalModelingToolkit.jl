# # Fabry-Perot AC transfer function
#
# Transfer functions through the ecosystem path: numeric linearization at the
# operating point (`linearize_proper`), then ControlSystemsBase (`ss`,
# `freqresp`, `poles`, `dcgain`). OMTK writes no frequency-sweep code of its
# own. The reference is the representation-parameterized analytic response
#
# ```math
# H(s) = \frac{t_1}{1 - r_1 r_2 P(s)^2}, \qquad P(s) = \mathrm{AP}_n(s)
# ```
#
# with ``\mathrm{AP}_n`` the order-``n`` all-pass cascade (here ``n = 1``), for a
# cavity on resonance (zero tuning).
#
# The input is the laser's quadrature modulation pair `(δA, δQ)`; the outputs
# are the quadratures of the intracavity field leaving `m1.p1`.

using Test
using OpticalModelingToolkit
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ControlSystemsBase: ss, freqresp, poles, dcgain

@testset "Fabry-Perot AC — ecosystem path vs analytic reference" begin
    ## Numerical parameters
    R1 = 0.99
    R2 = 0.999
    r1_v = sqrt(R1)
    r2_v = sqrt(R2)
    t1_v = sqrt(1 - R1)
    c_v = 299792458.0
    L_v = 0.1                       # resonant: zero tuning
    τ_v = L_v / c_v

    ## Model: PadeDelay(1), one all-pass section per direction, 4 states
    @named laser = LaserModulated(A = 1.0, ϕ = 0.0)
    @named m1 = Mirror(R = R1)
    @named space = FreeSpace(L = L_v, rep = PadeDelay(1))
    @named m2 = Mirror(R = R2)
    @named vac = Vacuum()
    eqs = [
        link(laser.p1, m1.p2)
        link(m1.p1, space.p1)
        link(space.p2, m2.p1)
        link(m2.p2, vac.p1)
    ]
    @named fp = System(eqs, t; systems = [laser, m1, space, m2, vac])

    inputs = [laser.δA, laser.δQ]
    outputs = [m1.p1.o.a[1], m1.p1.o.a[2]]

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
    ## P(s) = (1 - sτ/2)/(1 + sτ/2): all-pass, so |H| is exact at every
    ## frequency and only the phase carries the (ωτ)³/12 section error.
    Hana = function (ω)
        ap = (1 - im * ω * τ_v / 2) / (1 + im * ω * τ_v / 2)
        t1_v / (1 - r1_v * r2_v * ap^2)
    end
    for (k, ω) in enumerate(ωs)
        @test Hnum[k] ≈ Hana(ω) rtol = 1.0e-10
    end

    ## DC gain: the algebraic intracavity buildup t1/(1 - r1 r2), identical in
    ## every representation
    dc = dcgain(sys)
    @test dc[1, 1] ≈ t1_v / (1 - r1_v * r2_v) rtol = 1.0e-10
    @test Hnum[1] ≈ Hana(ωs[1]) rtol = 1.0e-10

    println(
        "✓ Fabry-Perot AC: ecosystem path matches the analytic reference over ",
        length(ωs), " points; |H(0)| = ",
        round(abs(dc[1, 1]), digits = 3)
    )
end
