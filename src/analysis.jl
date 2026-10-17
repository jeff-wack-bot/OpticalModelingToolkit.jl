# Linearization for the AC path (PLAN_rationalize_late §4).
#
# `ModelingToolkit.linearize` keeps torn algebraic unknowns (the circulating-
# field loops that survive `mtkcompile`) as states in the returned realization:
# it differentiates their constraints, which (a) demands input-derivative
# columns whenever the loop has direct feedthrough from an input — true for
# every all-pass PadeDelay network — and (b) leaves one spurious zero
# eigenvalue per algebraic unknown, so the realization is singular at s = 0
# and `dcgain` is wrong. For an index-1 DAE the honest reduction is static
# elimination of the algebraic block, which is what `linearize_proper` does
# with the Jacobians MTK already computes:
#
#     z = -g_z⁻¹ (g_x x + g_u u)   ⇒   A = f_x + f_z z_x,  B = f_u + f_z z_u,
#                                      C = h_x + h_z z_x,  D = h_u + h_z z_u.

"""
    linearize_proper(model, inputs, outputs; op = Dict(), t = 0.0, kwargs...)
        -> (; A, B, C, D), simplified_sys

Numerically linearize `model` between `inputs` and `outputs` at the operating
point `op`, returning a *proper* state-space realization whose states are
exactly the differential unknowns: index-1 algebraic variables (torn
circulating-field loops) are eliminated statically instead of being kept as
states behind differentiated constraints, so no input-derivative columns
appear and the realization is regular at `s = 0` (`dcgain` is exact).

Feed the matrices straight into `ControlSystemsBase`:

    mats, _ = linearize_proper(model, [laser.δA, laser.δQ], [pd.P])
    sys = ss(mats...)
    freqresp(sys, ωs)    # bodeplot, poles, dcgain, …

Keyword arguments are forwarded to `ModelingToolkit.linearization_function`.
"""
function linearize_proper(
        model::System, inputs, outputs;
        op = Dict(), t = 0.0, kwargs...
    )
    linfun, ssys = ModelingToolkit.linearization_function(
        model, inputs, outputs;
        op, warn_empty_op = false, kwargs...
    )
    prob = ModelingToolkit.LinearizationProblem(linfun, t)
    linres = prob.f(state_values(prob), parameter_values(prob), prob.t)
    (; f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u) = linres

    isempty(g_z) && return (A = f_x, B = f_u, C = h_x, D = h_u), ssys

    gz = lu(g_z; check = false)
    issuccess(gz) ||
        error(
        "g_z is not invertible: the DAE is of index > 1, or an " *
            "algebraic loop is singular at this operating point"
    )
    zx = -(gz \ g_x)
    zu = -(gz \ g_u)
    return (
            A = f_x + f_z * zx, B = f_u + f_z * zu,
            C = h_x + h_z * zx, D = h_u + h_z * zu,
        ), ssys
end
