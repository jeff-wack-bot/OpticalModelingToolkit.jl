# Electronics: linear filters for control loops.
#
# A servo in an interferometer is a linear filter between a sensor and an
# actuator. `ZPK` writes one as a state-space system, so it sits in the same MTK
# model as the optics and the suspension.

"""
    ZPK(; name, z = Float64[], p = Float64[], k = nothing, gain = 1.0)

A filter with transfer function

```math
H(s) = \\mathrm{gain}\\cdot k\\,\\frac{\\prod_i (s - z_i)}{\\prod_j (s - p_j)},
```

with zeros `z` and poles `p` in rad/s. Input `u`, output `y`.

When `k` is `nothing` it is ``\\prod_j |p_j| / \\prod_i |z_i|``, which gives unit
gain at DC for a filter with no zero at the origin, and `1` if a zero is at the
origin. This is Finesse's `zpk`, so a KatScript `zpk` line means the same filter
in both codes.

The roots must be real, and the filter proper (no more zeros than poles). It is
realized in controllable canonical form, one state per pole.
"""
@component function ZPK(; name, z = Float64[], p = Float64[], k = nothing, gain = 1.0)
    length(z) <= length(p) ||
        throw(ArgumentError("ZPK `$name`: $(length(z)) zeros and $(length(p)) poles is not proper"))
    if k === nothing
        pz = prod(abs, z; init = 1.0)
        k = pz == 0 ? 1.0 : prod(abs, p; init = 1.0) / pz
    end
    G = gain * k
    den = real_poly(p)                 # monic, ascending: den[end] == 1
    num = G .* real_poly(z)
    n = length(p)
    d = length(z) == n ? num[end] : 0.0
    rem = [i <= length(num) ? num[i] : 0.0 for i in 1:n] .- d .* den[1:n]

    @variables u(t) [input = true] y(t) [output = true]
    if n == 0
        return System([y ~ d * u], t, [u, y], []; name)
    end
    xs = [only(@variables $(Symbol(:x, i))(t) = 0.0) for i in 1:n]
    eqs = Equation[D(xs[i]) ~ xs[i + 1] for i in 1:(n - 1)]
    push!(eqs, D(xs[n]) ~ -sum(den[i] * xs[i] for i in 1:n) + u)
    push!(eqs, y ~ sum(rem[i] * xs[i] for i in 1:n) + d * u)
    return System(eqs, t, [xs; u; y], []; name)
end

"Real coefficients, lowest power first, of the monic polynomial with real `roots`."
function real_poly(roots)
    c = [1.0]
    for r in roots
        c = [0.0; c] .- r .* [c; 0.0]
    end
    return c
end
