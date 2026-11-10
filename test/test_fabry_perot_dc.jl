# # Fabry-Perot DC operating point
#
# The simplest meaningful network: a laser, two mirrors facing each other across
# a `FreeSpace`, and a photodetector on transmission. Every face-to-face contact
# is one `link`, which wires both directions. Each mirror's reflective face `p1`
# faces into the cavity.
#
# ```
#     Laser ── p2 [m1] p1 ══ space ══ p1 [m2] p2 ── pd
# ```
#
# With the algebraic propagation representation the whole system is algebraic,
# so `mtkcompile` produces an equation-only system that solves as a
# `NonlinearProblem`. The round-trip gain is ``r_1 r_2 e^{2i\Phi}``, where
# ``\Phi = \phi_1 + \phi_2`` is the total microscopic tuning, so the cavity is
# resonant at zero tuning and the field circulating out of `m1.p1` is
#
# ```math
# E_\mathrm{cav} = \frac{t_1 A}{1 - r_1 r_2 e^{2i\Phi}}.
# ```

using Test
using OpticalModelingToolkit
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using NonlinearSolve
using SciMLBase
using ControlSystemsBase: ss, dcgain

const R1, R2 = 0.99, 0.999
const r1, r2 = sqrt(R1), sqrt(R2)
const t1, t2 = sqrt(1 - R1), sqrt(1 - R2)
const λ0 = 1064.0e-9
const k0 = 2π / λ0

# A builder for the cavity. `m2` is a `MovableMirror`, so its position input
# must be bound: `x` pins it, and `x = nothing` leaves it free to be named as a
# linearization input.

function cavity(; x = 0.0, ϕ1 = 0.0, A = 1.0, rep = Algebraic())
    @named laser = Laser(; A)
    @named m1 = Mirror(R = R1, ϕ = ϕ1)
    @named space = FreeSpace(L = 0.1; rep)
    @named m2 = MovableMirror(R = R2, λ = λ0)
    @named pd = Photodetector()
    eqs = [
        link(laser.p1, m1.p2)
        link(m1.p1, space.p1)
        link(space.p2, m2.p1)
        link(m2.p2, pd.p1)
    ]
    x === nothing || push!(eqs, m2.x ~ x)
    @named fp = System(eqs, t; systems = [laser, m1, space, m2, pd])
    return fp, (; laser, m1, space, m2, pd)
end

# Zero is a valid starting point for this linear field network.

function solve_dc(model)
    sys = mtkcompile(model)
    prob = NonlinearProblem(sys, Dict(u => 0.0 for u in unknowns(sys)))
    sol = solve(prob, NewtonRaphson())
    @test SciMLBase.successful_retcode(sol)
    return sol
end

cfield(sol, port) = complex(sol[port.a[1]], sol[port.a[2]])

# ## Fields and energy conservation
#
# Detuned by ``\phi_1``, every field must match the analytic formulas, and
# reflected plus transmitted power must equal the input power.

@testset "fields and energy conservation" begin
    for ϕ1 in (0.0, 0.003, 0.4, π / 2)
        fp, c = cavity(; ϕ1)
        sol = solve_dc(fp)

        E_cav = t1 / (1 - r1 * r2 * cis(2ϕ1))          # out of m1.p1
        E_refl = -r1 * cis(-2ϕ1) + t1 * r2 * E_cav       # out of m1.p2
        E_trans = t2 * E_cav                             # out of m2.p2

        @test cfield(sol, c.m1.p1.o) ≈ E_cav rtol = 1.0e-12
        @test cfield(sol, c.m1.p2.o) ≈ E_refl rtol = 1.0e-12 atol = 1.0e-14
        @test cfield(sol, c.m2.p2.o) ≈ E_trans rtol = 1.0e-12
        @test sol[c.pd.P] ≈ abs2(E_trans) rtol = 1.0e-12
        @test abs2(cfield(sol, c.m1.p2.o)) + sol[c.pd.P] ≈ 1 rtol = 1.0e-12
    end
end

# ## The Airy fringe of a moving mirror
#
# Moving the end mirror by `x` tunes the round trip by ``2kx``. The transmitted
# power must trace
#
# ```math
# P_t(x) = \left| \frac{t_1 t_2 A}{1 - r_1 r_2 e^{2ikx}} \right|^2 .
# ```

P_t(x) = abs2(t1 * t2 / (1 - r1 * r2 * cis(2k0 * x)))

@testset "Airy fringe vs end-mirror position" begin
    for x in (0.0, 0.05e-9, -0.1e-9, 0.25e-9, 1.0e-9, λ0 / 4)
        fp, c = cavity(; x)
        @test solve_dc(fp)[c.pd.P] ≈ P_t(x) rtol = 1.0e-10
    end
end

# ## Radiation pressure
#
# `MovableMirror` reports the force of the light on both faces. On resonance the
# circulating power ``P_c`` peaks, and since nothing arrives from the detector
# side the end mirror feels ``F = (1 + R_2 - T_2) P_c / c = 2R_2 P_c/c``.

@testset "radiation pressure on the end mirror" begin
    c0 = 299792458.0
    fp, c = cavity(A = sqrt(2.0))       # 2 W input
    sol = solve_dc(fp)
    P_c = abs2(cfield(sol, c.m2.p1.i))
    @test P_c ≈ 2 * t1^2 / (1 - r1 * r2)^2 rtol = 1.0e-12
    @test sol[c.m2.F] ≈ 2 * R2 * P_c / c0 rtol = 1.0e-12
end

# ## The optical spring
#
# Off resonance the force depends on position, so the cavity acts as a spring
# with stiffness ``K = -\partial F/\partial x``. Linearizing the nonlinear model
# from `m2.x` to `m2.F` at a detuned operating point must reproduce the slope of
# the analytic force
#
# ```math
# F(x) = \frac{2 R_2}{c} \left| \frac{t_1 A}{1 - r_1 r_2 e^{2ikx}} \right|^2 .
# ```

@testset "optical spring from linearization" begin
    F(x) = 2R2 / 299792458.0 * abs2(t1 / (1 - r1 * r2 * cis(2k0 * x)))
    dF(x) = (h = 1.0e-14; (F(x + h) - F(x - h)) / 2h)
    fp, c = cavity(x = nothing)
    for x0 in (-2.0e-10, -5.0e-11, 0.0, 3.0e-11)
        mats, _ = linearize_proper(fp, [c.m2.x], [c.m2.F]; op = Dict(c.m2.x => x0))
        @test only(dcgain(ss(mats...))) ≈ dF(x0) rtol = 1.0e-5 atol = 1.0e-3
    end

    ## with delay states the linearization must still be taken at the steady
    ## state, where the delay lines hold the circulating field
    fp, c = cavity(x = nothing, rep = PadeDelay(2))
    mats, _ = linearize_proper(fp, [c.m2.x], [c.m2.F]; op = Dict(c.m2.x => -5.0e-11))
    @test only(dcgain(ss(mats...))) ≈ dF(-5.0e-11) rtol = 1.0e-5
end
