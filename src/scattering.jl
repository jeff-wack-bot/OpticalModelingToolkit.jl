# Lowering a composed optical model to the canonical transcendental scattering
# form (plans/PLAN_rationalize_late.md §5) — "rationalizing latest".
#
# For frequency-domain questions the model need not rationalize at all: every
# FreeSpace stays an exact delay `e^{-sτ}` with its static rotation folded into
# the adjacent matrix entries, exactly Etalon's `_collapse_linearspaces`. The
# result is Etalon's ScatteringSystem data (paper Eq. 11):
#
#   ( x'    )   ( M1  M2 ) ( E(s) x )
#   (       ) = (        ) (        )
#   ( y_out )   ( M3  M4 ) ( u_in   )
#
# with x indexing directed edges of the component graph, E(s) = diag(e^{-sτₖ}),
# and the exact transfer function G(s) = M3·(E(-s) − M1)⁻¹·M2 + M4.

"""
    to_scattering(model::System) -> NamedTuple

Lower a composed (un-`mtkcompile`d) optical model to the canonical
transcendental scattering form. Returns
`(M1, M2, M3, M4, delays, x, u, y)` — field-for-field the contents of an
Etalon `ScatteringSystem`; with Etalon loaded,
`Etalon.ScatteringSystem(model)` wraps this directly, giving exact
transcendental `freqresp` (FSR/comb physics included) and
`pade(sys; order)` as the late rational transformation.

The lowering works on the *scalar complex-amplitude* field: `x` indexes
directed edges of the component graph, each `FreeSpace` collapses into an
exact delay `τ = L/c` on its edge with the rotation `e^{i2πL/λ}` folded into
the consuming matrix entries. Any `PropagationRep` on a `FreeSpace` is
ignored — the delay here is exact by construction. Sources and sinks
(`Laser`, `Vacuum`, `Photodetector`) become the external channels `u`/`y`
(one input and one output per component, named after it); a `MovableMirror`
scatters at its reference position `x = 0`.

Requires a flat model: all components and `connect` equations at the top
level, every optical port wired.
"""
function to_scattering(model::System)
    comps = ModelingToolkit.get_systems(model)
    infoof = Dict{Symbol, NamedTuple}()
    for comp in comps
        info = SymbolicUtils.getmetadata(comp, OpticalComponent, nothing)
        info === nothing && error(
            "component $(nameof(comp)) carries no OpticalComponent metadata; " *
            "to_scattering only understands OpticalModelingToolkit components")
        infoof[nameof(comp)] = info
    end

    # Directed links (src, dst, delay, phase) with endpoints (component, port),
    # from the model's connect equations. OMTK connectors are directed, so
    # every connection is out-port → in-port.
    Port = Tuple{Symbol, Symbol}
    links = Tuple{Port, Port, Float64, Float64}[]
    for eq in ModelingToolkit.get_eqs(model)
        rhs = SymbolicUtils.unwrap_const(Symbolics.unwrap(eq.rhs))
        rhs isa ModelingToolkit.Connection || continue
        endpoints = rhs.systems
        length(endpoints) == 2 || error(
            "expected pairwise directed connections, got $eq")
        pair = map(endpoints) do sys
            parts = split(string(nameof(sys)), '₊')
            length(parts) == 2 || error(
                "connector $(nameof(sys)) is not component₊port; " *
                "to_scattering requires a flat model")
            (Symbol(parts[1]), Symbol(parts[2]))
        end
        isout = endswith.(String.(last.(pair)), "out")
        count(isout) == 1 || error("connection $eq is not out-port → in-port")
        src, dst = isout[1] ? pair : reverse(pair)
        push!(links, (src, dst, 0.0, 0.0))
    end

    # Collapse FreeSpace components: the links into fs.port_in and out of
    # fs.port_out fuse into one edge accumulating delay L/c and rotation 2πL/λ.
    spaces = Set(name for (name, info) in infoof if info.kind == :FreeSpace)
    while (i = findfirst(l -> l[2][1] in spaces, links)) !== nothing
        src, (fs, _), τ1, φ1 = links[i]
        j = findfirst(l -> l[1][1] == fs, links)
        j === nothing && error("FreeSpace $fs has no outgoing connection")
        _, dst, τ2, φ2 = links[j]
        info = infoof[fs]
        fused = (src, dst, τ1 + τ2 + info.L / info.c, φ1 + φ2 + 2π * info.L / info.λ)
        deleteat!(links, sort([i, j]))
        push!(links, fused)
    end
    edges = links

    K = length(edges)
    src_of = Dict{Port, Int}()
    dst_of = Dict{Port, Int}()
    for (k, (src, dst, _, _)) in enumerate(edges)
        haskey(src_of, src) && error("port $src sources two connections")
        haskey(dst_of, dst) && error("port $dst sinks two connections")
        src_of[src] = k
        dst_of[dst] = k
    end
    function edge_at(table, port)
        k = get(table, port, 0)
        k == 0 && error("optical port $(port[1]).$(port[2]) is not wired")
        return k
    end

    # External channels: one input (its emitted field) and one output (the
    # field arriving at it) per source/sink component, in declaration order.
    ext = [nameof(c) for c in comps
           if infoof[nameof(c)].kind in (:Laser, :Vacuum, :Photodetector)]
    chan = Dict(name => e for (e, name) in enumerate(ext))
    E = length(ext)

    M1 = zeros(ComplexF64, K, K)
    M2 = zeros(ComplexF64, K, E)
    M3 = zeros(ComplexF64, E, K)
    M4 = zeros(ComplexF64, E, E)
    delays = [e[3] for e in edges]
    rot = [cis(e[4]) for e in edges]    # folded FreeSpace rotation per edge

    for comp in comps
        name = nameof(comp)
        info = infoof[name]
        if info.kind in (:Laser, :Vacuum, :Photodetector)
            M2[edge_at(src_of, (name, :port_out)), chan[name]] = 1
            k = edge_at(dst_of, (name, :port_in))
            M3[chan[name], k] = rot[k]
        elseif info.kind in (:Mirror, :MovableMirror)
            S = ((:b_out, :f_in, info.t), (:b_out, :b_in, info.r),
                 (:f_out, :f_in, -info.r), (:f_out, :b_in, info.t))
            for (outport, inport, v) in S
                k_out = edge_at(src_of, (name, outport))
                k_in = edge_at(dst_of, (name, inport))
                M1[k_out, k_in] += v * rot[k_in]
            end
        elseif info.kind == :FreeSpace
            # already collapsed into the edges
        else
            error("to_scattering does not know component kind $(info.kind)")
        end
    end

    xnames = [Symbol(src[1], '_', src[2], '→', dst[1], '_', dst[2])
              for (src, dst, _, _) in edges]

    return (M1 = M1, M2 = M2, M3 = M3, M4 = M4, delays = delays,
        x = xnames, u = ext, y = copy(ext))
end
