# # Propagation representations
#
# Three validations from PLAN_rationalize_late §7:
#
# 1. **All-pass unit test**: every `PadeDelay` section has ``|H(i\omega)| = 1``
#    to machine precision; phase error scales ``\propto n^{-2}`` at fixed
#    ``\omega\tau``.
# 2. **FP convergence**: Fabry-Perot transmission with `PadeDelay(n)` converges
#    monotonically to the exact transcendental response; the first FSR peak
#    lands at ``c/2L`` once ``n \gtrsim \omega\tau``.
# 3. **Representation invariance in-band**: `Algebraic` and `PadeDelay(n)` agree
#    in the audio band for a 0.1 m cavity.
#
# All frequency responses go through the ecosystem AC path (numeric
# linearization via `linearize_proper`, sweeps via `ControlSystemsBase.freqresp`)
# — OMTK ships no frequency-sweep code.

using Test
using OpticalModelingToolkit
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ControlSystemsBase: ss, freqresp

const c_v = 299792458.0
const λ_v = 1064.0e-9

# The complex field transfer function from the modulation quadrature pair
# `(δA, δQ)` to an output quadrature pair: with the carrier at ``\phi = 0`` a
# complex drive ``\delta a = \delta A + i\,\delta Q`` produces the response
# column ``H\cdot(1, i)^\mathsf{T}`` combined as ``\mathrm{re} + i\,\mathrm{im}``.
# For the field-linear networks here both quadrature columns carry the same
# complex TF, so driving `δA` alone suffices: ``h(\omega) = H[\mathrm{re}, \delta A]
# + i\,H[\mathrm{im}, \delta A]``.

function field_tf(model, inputs, outputs, ωs)
    mats, _ = linearize_proper(model, inputs, outputs)
    resp = freqresp(ss(mats...), ωs)
    return [resp[1, 1, k] + im * resp[2, 1, k] for k in eachindex(ωs)]
end

# The exact all-pass section cascade, for realization checks:

ap_exact(ω, θ, n) = ((1 - im * ω * θ / 2) / (1 + im * ω * θ / 2))^n

# ## All-pass sections: unit magnitude, ``n^{-2}`` phase error
#
# A single `FreeSpace` between a modulated laser and a vacuum termination: the
# transfer function from the laser to the far side of the space is the delay
# realization itself. Each `PadeDelay` section is exactly all-pass, and the phase error at fixed
# ``\omega\tau`` shrinks by ``\times 4`` each time the order doubles.

@testset "all-pass sections: unit magnitude, n⁻² phase error" begin
    L_v = 0.1
    τ_v = L_v / c_v
    φ_v = 2π * L_v / λ_v

    orders = (1, 2, 4, 8, 16)
    ωτ_grid = (0.05, 0.2, 0.5, 1.0, 2.0, π, 5.0)
    phase_err = Dict{Int, Float64}()   # at the reference point ωτ = 0.5
    for n in orders
        @named laser = LaserModulated(A = 1.0, ϕ = 0.0)
        @named fs = FreeSpace(L = L_v, λ = λ_v, rep = PadeDelay(n))
        @named vac = Vacuum()
        eqs = [link(laser.p1, fs.p1); link(fs.p2, vac.p1)]
        @named line = System(eqs, t; systems = [laser, fs, vac])

        ωs = collect(ωτ_grid) ./ τ_v
        h = field_tf(
            line, [laser.δA, laser.δQ],
            [fs.p2.o.a[1], fs.p2.o.a[2]], ωs
        )

        for (k, ωτ) in enumerate(ωτ_grid)
            ## passivity: exactly all-pass at every order and frequency
            @test abs(h[k]) ≈ 1 atol = 1.0e-12
            ## the state-space realization is exactly the section cascade
            @test h[k] ≈ cis(φ_v) * ap_exact(ωτ / τ_v, τ_v / n, n) rtol = 1.0e-9
        end

        ω_ref = 0.5 / τ_v
        h_ref = field_tf(
            line, [laser.δA, laser.δQ],
            [fs.p2.o.a[1], fs.p2.o.a[2]], [ω_ref]
        )[1]
        phase_err[n] = abs(angle(h_ref * conj(cis(φ_v - ω_ref * τ_v))))
    end

    ## phase error ∝ n⁻² at fixed ωτ: each doubling of n gains ×4, and the
    ## absolute error matches the (ωτ)³/(12n²) estimate
    for n in (1, 2, 4, 8)
        @test phase_err[n] / phase_err[2n] ≈ 4 rtol = 0.15
    end
    @test phase_err[1] ≈ 0.5^3 / 12 rtol = 0.1
end

# ## Fabry-Perot convergence to the exact delay response
#
# On-resonance Fabry-Perot (``e^{2ikL} = -1``) with a modulated laser. The
# helper returns the composed model and the pieces needed for transfer
# functions.

function fp_model(L_v, rep)
    @named laser = LaserModulated(A = 1.0, ϕ = 0.0)
    @named m1 = Mirror(R = 0.99)
    @named space = FreeSpace(L = L_v, λ = λ_v, rep = rep)
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

# The exact transcendental transmission ``t_1 t_2 P / (1 + r_1 r_2 P^2)``,
# ``P = e^{i(\phi - \Omega\tau)}``, from the laser modulation to `m2.p2.o` (same
# scattering convention as the components, so ``e^{2i\phi} = -1`` on resonance):

function fp_exact_transmission(L_v)
    r1, r2 = sqrt(0.99), sqrt(0.999)
    t1, t2 = sqrt(1 - 0.99), sqrt(1 - 0.999)
    τ = L_v / c_v
    φ = 2π * L_v / λ_v
    return ω -> begin
        P = cis(φ - ω * τ)
        t1 * t2 * P / (1 + r1 * r2 * P^2)
    end
end

# Resonant macroscopic length ≈ 0.1 m: ``L = (2m+1)\lambda/4 \Rightarrow
# e^{2ikL} = -1``.

const L_res = (2 * round(Int, 0.1 / (λ_v / 4)) + 1) * λ_v / 4

@testset "Fabry-Perot convergence to the exact delay response" begin
    τ_v = L_res / c_v
    fsr = π / τ_v                    # FSR in angular frequency, = 2π·c/2L
    ## Band up to 60% of the FSR: away from the resonance peaks the response
    ## is smooth and the pointwise relative error is dominated by the phase
    ## error of the representation. (On the peaks themselves, any peak shift
    ## gives O(1) pointwise error at every order — peak location is tested
    ## separately below.)
    ωs = collect(range(0.01 * fsr, 0.6 * fsr; length = 121))
    H_ex = fp_exact_transmission(L_res)

    maxerr = Float64[]
    orders = (1, 2, 4, 8, 16)
    for n in orders
        fp, laser, m1, m2 = fp_model(L_res, PadeDelay(n))
        h = field_tf(
            fp, [laser.δA, laser.δQ],
            [m2.p2.o.a[1], m2.p2.o.a[2]], ωs
        )
        push!(maxerr, maximum(abs.(h .- H_ex.(ωs)) ./ abs.(H_ex.(ωs))))
    end
    ## the approximation order is a convergence knob: error decreases
    ## monotonically with n over the fixed band
    @test all(diff(maxerr) .< 0)
    @test maxerr[end] < 0.01

    ## FSR physics: the first free-spectral-range peak sits at Ω = π/τ = 2π·c/2L
    fp, laser, m1, m2 = fp_model(L_res, PadeDelay(16))
    ω_scan = collect(range(0.9 * fsr, 1.1 * fsr; length = 401))
    h_scan = field_tf(
        fp, [laser.δA, laser.δQ],
        [m2.p2.o.a[1], m2.p2.o.a[2]], ω_scan
    )
    ω_peak = ω_scan[argmax(abs.(h_scan))]
    @test abs(ω_peak - fsr) / fsr < 0.01

    println(
        "✓ FP convergence: max rel err by order ",
        Dict(zip(orders, round.(maxerr, sigdigits = 2))),
        "  FSR peak off by ", round(100 * abs(ω_peak - fsr) / fsr, sigdigits = 2), "%"
    )
end

# ## Representation invariance in the audio band
#
# `PadeDelay` is audio-band-exact at any order (phase error ``(\omega\tau)^3/12n^2
# \lesssim 10^{-13}`` rad here); the adiabatic `Algebraic` representation agrees
# to the cavity-pole correction ``\sim 2 r_1 r_2 \omega\tau/(1 - r_1 r_2)``,
# ``\approx 1\%`` at 10 kHz for this 0.1 m cavity.

@testset "representation invariance in the audio band" begin
    τ_v = L_res / c_v
    H_ex = fp_exact_transmission(L_res)
    ω_audio = 2π .* (10.0 .^ range(1, 4; length = 13))   # 10 Hz – 10 kHz

    fp_alg, laser, m1, m2 = fp_model(L_res, Algebraic())
    h_alg = field_tf(
        fp_alg, [laser.δA, laser.δQ],
        [m2.p2.o.a[1], m2.p2.o.a[2]], ω_audio
    )

    fp_p1, laser, m1, m2 = fp_model(L_res, PadeDelay(1))
    h_p1 = field_tf(
        fp_p1, [laser.δA, laser.δQ],
        [m2.p2.o.a[1], m2.p2.o.a[2]], ω_audio
    )

    fp_p4, laser, m1, m2 = fp_model(L_res, PadeDelay(4))
    h_p4 = field_tf(
        fp_p4, [laser.δA, laser.δQ],
        [m2.p2.o.a[1], m2.p2.o.a[2]], ω_audio
    )

    for (k, ω) in enumerate(ω_audio)
        ## PadeDelay is audio-band-exact at any order; the tolerance is the
        ## finesse-amplified roundoff floor of the closed loop, not a
        ## representation error
        @test h_p1[k] ≈ H_ex(ω) rtol = 1.0e-6
        @test h_p4[k] ≈ H_ex(ω) rtol = 1.0e-6
        ## the adiabatic representation agrees to the cavity-pole correction
        @test h_alg[k] ≈ H_ex(ω) rtol = 0.02
    end
    ## ... and tightly at the bottom of the band
    @test h_alg[1] ≈ H_ex(ω_audio[1]) rtol = 1.0e-4
end
