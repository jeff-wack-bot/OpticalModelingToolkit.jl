# # The transcendental scattering bridge
#
# `to_scattering` lowers the composed model to the canonical scattering form
# with every `FreeSpace` an *exact* delay (PLAN_rationalize_late §5, validation
# 5). Three layers of validation:
#
# 1. the exact scalar transfer function reproduces the analytic Fabry-Perot
#    response, transcendentally in ``\Omega`` (FSR physics included);
# 2. the **A-matrix keystone test**: OMTK's linearized two-photon (quadrature)
#    response equals ``A\cdot\mathrm{diag}(g(i\Omega), \overline{g(-i\Omega)})\cdot
#    A^{-1}`` blockwise per port pair, on an undetuned and a detuned cavity,
#    where ``g`` is the scalar network evaluated with the *matching* delay
#    representation (the identity is exact for any real-rational representation)
#    — this cross-validates OMTK's realification, its linearization path, and
#    `to_scattering` against each other;
# 3. **rep convergence across the bridge**: a high-order `PadeDelay` model tracks
#    the exact-delay scattering response within the known phase-error bound.
#
# When [Etalon](https://github.com) is available (e.g. dev'd into the
# environment), the same matrices drive Etalon's `ScatteringSystem` solver via
# the package extension, closing the cross-tool validation web.

using Test
using OpticalModelingToolkit
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ControlSystemsBase: ss, freqresp
using LinearAlgebra

const λ_sc = 1064.0e-9
const c_sc = 299792458.0

# Transfer matrix ``G(s) = M_3 (E(-s) - M_1)^{-1} M_2 + M_4`` of a
# `to_scattering` result. `delayfun(s, τ)` evaluates one edge's propagation
# dynamics: `exp(-s*τ)` is the exact transcendental form (what Etalon's
# `ScatteringSystem` evaluates); an all-pass rational stands in for a `PadeDelay`
# time-domain representation.

function scattering_tf(parts, s; delayfun = (s, τ) -> exp(-s * τ))
    Einv = Diagonal([1 / delayfun(s, τ) for τ in parts.delays])
    return parts.M3 * ((Einv - parts.M1) \ parts.M2) + parts.M4
end

# The order-`n` all-pass cascade realized by `PadeDelay(n)`:

ap_delay(n) = (s, τ) -> ((1 - s * τ / 2n) / (1 + s * τ / 2n))^n

# Fabry-Perot with a modulated laser; `L` sets the (de)tuning.

function fp_model(L_v, rep = Algebraic())
    @named laser = LaserModulated(A = 1.0, ϕ = 0.0)
    @named m1 = Mirror(R = 0.99)
    @named space = FreeSpace(L = L_v, λ = λ_sc, rep = rep)
    @named m2 = Mirror(R = 0.999)
    @named vac = Vacuum()
    eqs = [
        link(laser.p1, m1.p1)
        link(m1.p2, space.p1)
        link(space.p2, m2.p1)
        link(m2.p2, vac.p1)
    ]
    @named fp = System(eqs, t; systems = [laser, m1, space, m2, vac])
    return fp, laser, m1, m2
end

# Resonant macroscopic length ≈ 0.1 m: ``e^{2ikL} = -1``.

const L_res_sc = (2 * round(Int, 0.1 / (λ_sc / 4)) + 1) * λ_sc / 4

# ## Exact transcendental Fabry-Perot response
#
# `to_scattering` must expose six directed edges and preserve the FreeSpace
# delay ``\tau = L/c`` exactly on the two intracavity edges. The resulting
# scalar transfer function reproduces the analytic transmission and reflection,
# swept transcendentally across two FSRs with no rationalization anywhere.

@testset "to_scattering: exact transcendental Fabry-Perot response" begin
    fp, laser, m1, m2 = fp_model(L_res_sc)
    parts = to_scattering(fp)

    ## 6 directed edges: laser↔m1 (2), m1↔m2 through the spaces (2), m2↔vac (2)
    @test length(parts.delays) == 6
    @test length(parts.x) == 6
    @test sort(parts.u) == sort(parts.y)
    ## the FreeSpace delay τ = L/c survives exactly on the two intracavity edges
    τ_v = L_res_sc / c_sc
    @test count(≈(τ_v), parts.delays) == 2
    @test count(==(0.0), parts.delays) == 4

    iu = findfirst(==(:laser), parts.u)
    iy_t = findfirst(==(:vac), parts.y)     # transmission arrives at the vacuum
    iy_r = findfirst(==(:laser), parts.y)   # reflection arrives back at the laser

    r1, r2 = sqrt(0.99), sqrt(0.999)
    t1, t2 = sqrt(1 - 0.99), sqrt(1 - 0.999)
    φ = 2π * L_res_sc / λ_sc

    ## Transcendental sweep across two FSRs. Tolerance: near resonance the
    ## closed-loop denominator amplifies the large-angle trig roundoff of
    ## φ ≈ 6e5 rad by the finesse (~2e2), so ~1e-8 relative is the float-level
    ## agreement, not a physics error.
    for f in 10.0 .^ range(2, log10(2.2 * c_sc / (2 * L_res_sc)); length = 40)
        Ω = 2π * f
        P = cis(φ - Ω * τ_v)
        G = scattering_tf(parts, im * Ω)
        @test G[iy_t, iu] ≈ t1 * t2 * P / (1 + r1 * r2 * P^2) rtol = 1.0e-6
        @test G[iy_r, iu] ≈ r1 + t1^2 * r2 * P^2 / (1 + r1 * r2 * P^2) rtol = 1.0e-6
    end
end

# ## A-matrix keystone: two-photon block ↔ scalar sidebands
#
# The realification ``(\mathrm{re}, \mathrm{im}) = T (a, \bar a)``,
# ``T = \tfrac12 \begin{bmatrix} 1 & 1 \\ -i & i \end{bmatrix}``, ties OMTK's
# linearized quadrature blocks to the rep-matched scalar network. Because the
# representation is a real-rational (conj-symmetric) function of ``s``, the
# blockwise identity holds exactly at every frequency — checked here on both an
# undetuned and a detuned cavity, in both the `Algebraic` and `PadeDelay(2)`
# representations.

@testset "A-matrix keystone: two-photon block ↔ scalar sidebands" begin
    T = [1 1; -im im] / 2

    for (label, L_v) in (("undetuned", L_res_sc), ("detuned", L_res_sc + λ_sc / 17)),
            (rep, delayfun) in ((Algebraic(), (s, τ) -> 1.0), (PadeDelay(2), ap_delay(2)))

        fp, laser, m1, m2 = fp_model(L_v, rep)
        parts = to_scattering(fp)
        iu = findfirst(==(:laser), parts.u)

        ## OMTK linearized quadrature blocks vs the rep-matched scalar network:
        ## exact blockwise identity at every frequency, since the rep is a
        ## real-rational function of s (conj-symmetric)
        Ωs = 2π .* (10.0 .^ range(3, 9.3; length = 9))
        for (out, ychan) in ((m2.p2.o, :vac), (m1.p1.o, :laser))
            mats, _ = linearize_proper(fp, [laser.δA, laser.δQ], [out.a[1], out.a[2]])
            resp = freqresp(ss(mats...), Ωs)
            iy = findfirst(==(ychan), parts.y)
            for (k, Ω) in enumerate(Ωs)
                g₊ = scattering_tf(parts, im * Ω; delayfun)[iy, iu]
                g₋ = conj(scattering_tf(parts, -im * Ω; delayfun)[iy, iu])
                pred = T * Diagonal([g₊, g₋]) / T
                @test resp[:, :, k] ≈ pred rtol = 1.0e-6 atol = 1.0e-12
            end
        end
    end
    println("✓ A-matrix keystone holds, undetuned and detuned, both reps")
end

# ## Rep convergence across the bridge
#
# A `PadeDelay(16)` time-domain model tracks the exact-delay scattering response
# within the ``(\omega\tau)^3/12n^2`` phase-error bound out to half an FSR.

@testset "rep convergence across the bridge" begin
    fp, laser, m1, m2 = fp_model(L_res_sc, PadeDelay(16))
    parts = to_scattering(fp)
    iu = findfirst(==(:laser), parts.u)
    iy = findfirst(==(:vac), parts.y)

    Ωs = 2π .* (10.0 .^ range(5, log10(0.5 * c_sc / (2 * L_res_sc)); length = 12))
    mats, _ = linearize_proper(fp, [laser.δA, laser.δQ], [m2.p2.o.a[1], m2.p2.o.a[2]])
    resp = freqresp(ss(mats...), Ωs)
    for (k, Ω) in enumerate(Ωs)
        h_omtk = resp[1, 1, k] + im * resp[2, 1, k]
        h_exact = scattering_tf(parts, im * Ω)[iy, iu]
        @test h_omtk ≈ h_exact rtol = 2.0e-3
    end
end

# ## Etalon cross-check through the extension
#
# The cross-tool leg: when Etalon is loadable, its `ScatteringSystem` must
# evaluate the identical transcendental response through the package extension.

if Base.find_package("Etalon") !== nothing
    using Etalon: Etalon
    @testset "Etalon cross-check through the extension" begin
        fp, laser, m1, m2 = fp_model(L_res_sc)
        parts = to_scattering(fp)
        sys = Etalon.ScatteringSystem(fp)
        @test sys isa Etalon.ScatteringSystem
        for f in (1.0e3, 1.0e6, 1.0e9, 2.3e9)
            s = im * 2π * f
            @test sys(s) ≈ scattering_tf(parts, s) rtol = 1.0e-12
        end
        println("✓ Etalon ScatteringSystem matches to_scattering exactly")
    end
else
    @info "Etalon not in the test environment — cross-tool leg skipped " *
        "(dev ~/dev/Etalon into the env to enable it)"
end
