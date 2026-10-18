# # Optomechanics: optical spring and lock acquisition
#
# A Fabry-Perot cavity whose end mirror hangs on a suspension, so optics,
# mechanics and control share one MTK model. The suspension is written directly
# as equations: a mass `M` on a spring of resonance `fm` and quality factor
# `Q`, pushed by the mirror's radiation pressure force `F` and an external
# force `Fext`.
#
# ```
#     Laser ── m1 ══ space ══ m2 ── pd
#                              ⟂
#                         suspension (M, fm, Q)
# ```

using Test
using OpticalModelingToolkit
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ControlSystemsBase: ss, poles
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using SciMLBase

const λ_om = 1064.0e-9
const c_om = 299792458.0
const R1_om, R2_om = 0.99, 0.999
const M_om, fm_om, Q_om = 0.1, 20.0, 2.0
const k_om = M_om * (2π * fm_om)^2

# The radiation pressure force on the end mirror as a function of its position
# `x` (zero on resonance), for input power `Pin`:

function F_analytic(x, Pin)
    r1, r2, t1 = sqrt(R1_om), sqrt(R2_om), sqrt(1 - R1_om)
    return 2R2_om / c_om * Pin * abs2(t1 / (1 - r1 * r2 * cis(2 * (2π / λ_om) * x)))
end

# The model. `controller` adds equations closing a loop around the plant; with
# none, `Fext` is a free input.

function suspended_cavity(; Pin, x_rest, controller = nothing)
    @named laser = Laser(A = sqrt(Pin))
    @named m1 = Mirror(R = R1_om)
    @named space = FreeSpace(L = 0.1)
    @named m2 = MovableMirror(R = R2_om, λ = λ_om)
    @named pd = Photodetector()
    @variables x(t) = x_rest v(t) = 0.0 Fext(t) [input = true]
    ωm = 2π * fm_om
    eqs = [
        link(laser.p1, m1.p2)
        link(m1.p1, space.p1)
        link(space.p2, m2.p1)
        link(m2.p2, pd.p1)
        m2.x ~ x
        D(x) ~ v
        M_om * D(v) ~ -k_om * (x - x_rest) - (M_om * ωm / Q_om) * v + m2.F + Fext
    ]
    controller === nothing || append!(eqs, controller(pd, Fext))
    @named plant = System(eqs, t; systems = [laser, m1, space, m2, pd])
    return plant, (; x, v, Fext, m2, pd)
end

# ## The optical spring
#
# Detuned from resonance, the force depends on position, so radiation pressure
# acts as a spring of stiffness ``K = -\partial F/\partial x`` in parallel with
# the suspension. The linearized opto-mechanical model, from the external force
# to the position, must have exactly the poles of
#
# ```math
# M s^2 + \frac{M\omega_m}{Q} s + (k + K).
# ```
#
# The equilibrium position balances the suspension and the light,
# ``k(x^* - x_\mathrm{rest}) = F(x^*)``. Beyond the resonance (``x^* > 0``)
# the light stiffens the suspension. Short of it the optical spring is negative
# and much stronger than the suspension, so the mirror is unstable.

@testset "optical spring shifts the suspension resonance" begin
    Pin = 50.0
    for x_star in (0.3e-9, -0.3e-9)
        x_rest = x_star - F_analytic(x_star, Pin) / k_om
        plant, c = suspended_cavity(; Pin, x_rest)
        mats, _ = linearize_proper(plant, [c.Fext], [c.x]; op = Dict(c.x => x_star, c.v => 0.0))
        h = 1.0e-14
        K = -(F_analytic(x_star + h, Pin) - F_analytic(x_star - h, Pin)) / 2h
        ## analytic poles of M s² + (Mωm/Q) s + (k + K)
        a, b, cc = M_om, M_om * 2π * fm_om / Q_om, k_om + K
        disc = sqrt(complex(b^2 - 4a * cc))
        expected = [(-b + disc) / 2a, (-b - disc) / 2a]
        p = poles(ss(mats...))
        @test sort(p; by = z -> (imag(z), real(z))) ≈
            sort(expected; by = z -> (imag(z), real(z))) rtol = 1.0e-6
        if x_star > 0
            ## stiffened: the 20 Hz resonance moves up by more than 2×
            @test sqrt(cc / M_om) / 2π > 2fm_om
        else
            ## anti-spring: unstable
            @test maximum(real, p) > 0
        end
    end
end

# ## Lock acquisition in the time domain
#
# Now close a loop: an integrating servo pushes the mirror with ``F_\mathrm{ext}
# = g_I\int (P_\mathrm{set} - P)\,dt`` to hold the transmitted power at half its
# peak. The mirror starts at rest half a linewidth away from the lock point, on
# the low-power side; the servo pulls it in and holds it there. The model is a
# mass-matrix DAE (the field network stays algebraic), so it needs a Rosenbrock
# solver.

@testset "lock acquisition" begin
    Pin = 1.0
    r1r2 = sqrt(R1_om * R2_om)
    P_peak = (1 - R1_om) * (1 - R2_om) / (1 - r1r2)^2 * Pin
    P_set = P_peak / 2
    ## half-power point on the far side of the resonance, where dP/dx < 0
    P_t(x) = Pin * (1 - R1_om) * (1 - R2_om) /
        abs2(1 - r1r2 * cis(2 * (2π / λ_om) * x))
    lo, hi = 0.0, λ_om / 8
    for _ in 1:100
        mid = (lo + hi) / 2
        P_t(mid) > P_set ? (lo = mid) : (hi = mid)
    end
    x_lock = lo
    slope = (P_t(x_lock + 1.0e-13) - P_t(x_lock - 1.0e-13)) / 2.0e-13
    gI = 2π * 5.0 * k_om / abs(slope)            # ~5 Hz unity-gain frequency

    function servo(pd, Fext)
        @variables u(t) = 0.0
        return [D(u) ~ -gI * (P_set - pd.P), Fext ~ u]
    end
    linewidth = λ_om / (4π) * (1 - r1r2) / sqrt(r1r2)
    plant, c = suspended_cavity(; Pin, x_rest = x_lock + 0.5 * linewidth, controller = servo)
    sys = mtkcompile(plant)
    prob = ODEProblem(sys, [], (0.0, 3.0))
    sol = solve(prob, Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-12)
    @test SciMLBase.successful_retcode(sol)
    ## starts well below the setpoint ...
    @test sol[c.pd.P][1] < 0.9P_set
    ## ... and ends locked on it
    @test sol(3.0; idxs = c.pd.P) ≈ P_set rtol = 1.0e-3
    @test sol(3.0; idxs = c.x) ≈ x_lock rtol = 1.0e-6
end
