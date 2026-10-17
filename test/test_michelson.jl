# # Michelson interferometer
#
# A `Beamsplitter` splits the laser into two arms closed by perfectly reflecting
# end mirrors; the returning beams recombine at the beamsplitter, and a
# photodetector watches the antisymmetric port.
#
# ```
#                  ey
#                  ║ Ly
#                  p2
#     Laser ── p1 [bs] p3 ══ Lx ══ ex
#                  p4
#                  │
#                  pd
# ```
#
# Light entering `p1` reflects out of `p2` (the y arm) and transmits out of `p3`
# (the x arm). Back at the beamsplitter, the y return transmits to `p4` and the x
# return reflects there with ``-r``, so the antisymmetric field is
#
# ```math
# E_\mathrm{AS} = r t A \left(e^{2ikL_y} - e^{2ikL_x}\right),
# \qquad P_\mathrm{AS} = 4 R T A^2 \sin^2 k(L_y - L_x).
# ```

using Test
using OpticalModelingToolkit
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using NonlinearSolve
using LinearAlgebra

const λ_mi = 1064.0e-9
const c_mi = 299792458.0

function michelson(; Lx, Ly, R = 0.5, laser = LaserModulated)
    @named src = laser()
    @named bs = Beamsplitter(; R)
    @named sx = FreeSpace(L = Lx, λ = λ_mi)
    @named sy = FreeSpace(L = Ly, λ = λ_mi)
    @named ex = Mirror(R = 1.0)
    @named ey = Mirror(R = 1.0)
    @named vx = Vacuum()
    @named vy = Vacuum()
    @named pd = Photodetector()
    eqs = [
        link(src.p1, bs.p1)
        link(bs.p2, sy.p1)
        link(sy.p2, ey.p1)
        link(ey.p2, vy.p1)
        link(bs.p3, sx.p1)
        link(sx.p2, ex.p1)
        link(ex.p2, vx.p1)
        link(bs.p4, pd.p1)
    ]
    @named mich = System(eqs, t; systems = [src, bs, sx, sy, ex, ey, vx, vy, pd])
    return mich, (; src, bs, pd)
end

# ## The DC fringe
#
# Stepping the arm-length difference through half a wavelength sweeps the
# antisymmetric port from dark to bright and back; the symmetric port (back
# toward the laser) carries the rest of the power.

@testset "Michelson DC fringe and energy conservation" begin
    Lx = 1.0
    for R in (0.5, 0.3), ΔL in range(0, λ_mi / 2; length = 7)
        mich, c = michelson(; Lx, Ly = Lx + ΔL, R, laser = Laser)
        sys = mtkcompile(mich)
        sol = solve(NonlinearProblem(sys, Dict(u => 0.0 for u in unknowns(sys))))
        k = 2π / λ_mi
        P_as = abs2(sqrt(R * (1 - R)) * (cis(2k * (Lx + ΔL)) - cis(2k * Lx)))
        ## arm phases 2kL ~ 1e7 rad: trig roundoff limits agreement to ~1e-9
        @test sol[c.pd.P] ≈ P_as atol = 1.0e-8
        P_sym = sol[c.src.p1.i.a[1]]^2 + sol[c.src.p1.i.a[2]]^2
        @test P_sym + sol[c.pd.P] ≈ 1 rtol = 1.0e-10
    end
end

# ## The exact antisymmetric-port response
#
# `to_scattering` needs no beamsplitter-specific code: it reads the scattering
# table each component stamps into its metadata. The exact transfer function
# from the laser to the detector, with arm delays ``\tau_{x,y} = L_{x,y}/c``, is
#
# ```math
# G(i\Omega) = r t \left(e^{2i(kL_y - \Omega\tau_y)} - e^{2i(kL_x - \Omega\tau_x)}\right).
# ```

@testset "to_scattering through a beamsplitter" begin
    Lx, Ly = 1.0, 1.3
    mich, _ = michelson(; Lx, Ly, R = 0.4)
    parts = to_scattering(mich)
    iu = findfirst(==(:src), parts.u)
    iy = findfirst(==(:pd), parts.y)
    k = 2π / λ_mi
    r, τ = sqrt(0.4), sqrt(0.6)
    for f in (0.0, 1.0e3, 1.0e6, 7.7e7, 3.0e8)
        s = im * 2π * f
        G = parts.M3 * ((Diagonal(exp.(s .* parts.delays)) - parts.M1) \ parts.M2) + parts.M4
        expected = r * τ * (exp(2im * k * Ly - 2s * Ly / c_mi) - exp(2im * k * Lx - 2s * Lx / c_mi))
        @test G[iy, iu] ≈ expected rtol = 1.0e-9
    end
end
