# Linearization for transfer functions.
#
# `ModelingToolkit.linearize` keeps algebraic unknowns (the circulating-field
# loops) as states in the returned realization: it differentiates their
# constraints, which (a) demands input-derivative columns whenever the loop has
# direct feedthrough from an input, true for every all-pass PadeDelay network,
# and (b) leaves one spurious zero eigenvalue per algebraic unknown, so the
# realization is singular at s = 0 and `dcgain` is wrong. For an index-1 DAE
# the honest reduction is static elimination of the algebraic block, which is
# what `linearize_proper` does with the Jacobians MTK already computes:
#
#     z = -g_z⁻¹ (g_x x + g_u u)   ⇒   A = f_x + f_z z_x,  B = f_u + f_z z_u,
#                                      C = h_x + h_z z_x,  D = h_u + h_z z_u.
#
# MTK's linearization solves only the algebraic unknowns at the operating
# point; differential states keep whatever `op` or their guesses say. For a
# `PadeDelay` network that means linearizing around a cavity whose delay lines
# hold zero field, which is wrong as soon as anything nonlinear (a
# `MovableMirror`, a photodetector) depends on the field. So the differential
# states that `op` does not pin are first driven to steady state by Newton's
# method on the reduced dynamics `ẋ = f(x, z(x))`, whose Jacobian is `A`.

"""
    linearize_proper(model, inputs, outputs; op = Dict(), t = 0.0, kwargs...)
        -> (; A, B, C, D), simplified_sys

Numerically linearize `model` between `inputs` and `outputs` at a steady state,
returning a *proper* state-space realization whose states are exactly the
differential unknowns. Index-1 algebraic variables (the circulating fields) are
eliminated statically instead of being kept as states behind differentiated
constraints, so no input-derivative columns appear and the realization is
regular at `s = 0` (`dcgain` is exact).

`op` sets inputs, parameters and any differential states to hold fixed; every
other differential state is solved for steady state (`ẋ = 0`). If those states
have no unique steady state (a free integrator, a free mass), pin them in `op`.

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
    build(op) = ModelingToolkit.linearization_function(
        model, inputs, outputs;
        op, warn_empty_op = false, warn_initialize_determined = false, kwargs...
    )
    linfun, ssys = build(op)
    pinned = Set(Symbolics.unwrap(k) for k in keys(op))
    free_states = filter(x -> Symbolics.unwrap(x) ∉ pinned, unknowns(ssys)[linfun.diff_idxs])
    if !isempty(free_states)
        # MTK leaves differential states without an initial value
        # unconstrained (least squares on their guesses), so rebuild with the
        # free states in `op`: they become `Initial` parameters that
        # `set_state!` can move.
        linfun, ssys = build(merge(Dict{Any, Any}(op), Dict(x => 0.0 for x in free_states)))
    end
    free = [i for i in linfun.diff_idxs if any(isequal(unknowns(ssys)[i]), free_states)]
    prob = ModelingToolkit.LinearizationProblem(linfun, t)
    evaluate() = prob.f(state_values(prob), parameter_values(prob), t)
    linres = evaluate()
    for _ in 1:10
        isempty(free) && break
        du = similar(linres.x)
        linfun.prob.f(du, linres.x, linres.p, t)
        A = eliminate_algebraic(linres).A
        rows = [findfirst(==(i), linfun.diff_idxs) for i in free]
        Afree = lu(A[rows, rows]; check = false)
        issuccess(Afree) || error(
            "no unique steady state for the differential states " *
                "$(unknowns(ssys)[free]); pin them in `op`"
        )
        step = Afree \ du[free]
        for (i, Δ) in zip(free, step)
            SymbolicIndexingInterface.set_state!(linfun, linres.x[i] - Δ, i)
        end
        linres = evaluate()
        maximum(abs, step) ≤ 1.0e-12 * max(1, maximum(abs, linres.x[free])) && break
    end
    return eliminate_algebraic(linres), ssys
end

function eliminate_algebraic(linres)
    (; f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u) = linres
    isempty(g_z) && return (A = f_x, B = f_u, C = h_x, D = h_u)
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
    )
end
