# Lowering a composed optical model to the canonical transcendental scattering
# form: "rationalizing latest".
#
# For frequency-domain questions the model need not rationalize at all: every
# FreeSpace stays an exact delay `e^{-sτ}` with its static rotation folded into
# the adjacent matrix entries (as Etalon's `_collapse_linearspaces` does). The
# result is the data of Etalon's `ScatteringSystem`:
#
#   ( x'    )   ( M1  M2 ) ( E(s) x )
#   (       ) = (        ) (        )
#   ( y_out )   ( M3  M4 ) ( u_in   )
#
# with x indexing directed edges of the component graph, E(s) = diag(e^{-sτₖ}),
# and the exact transfer function G(s) = M3·(E(-s) − M1)⁻¹·M2 + M4.

"""
    to_scattering(model::System) -> NamedTuple

Lower a composed (not yet `mtkcompile`d) optical model to the canonical
transcendental scattering form. Returns `(M1, M2, M3, M4, delays, x, u, y)`,
field for field the contents of an Etalon `ScatteringSystem`, so that

    G(s) = M3 * ((Diagonal(exp.(s .* delays)) - M1) \\ M2) + M4

is the exact transfer matrix from the external inputs `u` to the external
outputs `y`, FSR physics included. With Etalon loaded,
`Etalon.ScatteringSystem(model)` wraps this directly.

The lowering works on the scalar complex amplitude and reads only the
[`OpticalComponent`](@ref) metadata of each component: `x` indexes directed
edges, each [`Propagator`](@ref) (`FreeSpace`) collapses into an exact delay on
its edges with its rotation folded into the next matrix entry, and each
[`Terminal`](@ref) (`Laser`, `Vacuum`, `Photodetector`, …) contributes one
input and one output channel named after it. Any `PropagationRep` is ignored;
the delays here are exact. Components scatter at their constructor values (a
`MovableMirror` at `x = 0`).

Requires a flat model: all components and `link`/`connect` equations at the top
level, and every optical port wired.
"""
function to_scattering(model::System)
    comps = ModelingToolkit.get_systems(model)
    role = Dict{Symbol, Any}()
    for comp in comps
        r = SymbolicUtils.getmetadata(comp, OpticalComponent, nothing)
        r === nothing && error(
            "component $(nameof(comp)) carries no OpticalComponent metadata; " *
                "to_scattering only understands optical components"
        )
        role[nameof(comp)] = r
    end

    # Directed links (src, dst, delay, phase) between ports (component, port),
    # from the model's connect equations, each `comp₊port₊o` → `comp₊port₊i`.
    Port = Tuple{Symbol, Symbol}
    links = Tuple{Port, Port, Float64, Float64}[]
    for eq in ModelingToolkit.get_eqs(model)
        rhs = SymbolicUtils.unwrap_const(Symbolics.unwrap(eq.rhs))
        rhs isa ModelingToolkit.Connection || continue
        ends = map(rhs.systems) do sys
            parts = split(string(nameof(sys)), '₊')
            length(parts) == 3 || error(
                "connector $(nameof(sys)) is not component₊port₊(i|o); " *
                    "to_scattering requires a flat model of optical ports"
            )
            (Symbol(parts[1]), Symbol(parts[2]), parts[3])
        end
        length(ends) == 2 && sort([e[3] for e in ends]) == ["i", "o"] ||
            error("expected a pairwise o → i connection, got $eq")
        src, dst = ends[1][3] == "o" ? ends : reverse(ends)
        push!(links, (src[1:2], dst[1:2], 0.0, 0.0))
    end

    # Collapse propagators: the link into one side and the link out of the
    # other fuse into a single edge carrying the delay and rotation.
    for (name, r) in role
        r isa Propagator || continue
        for (a, b) in ((:p1, :p2), (:p2, :p1))
            i = findfirst(l -> l[2] == (name, a), links)
            j = findfirst(l -> l[1] == (name, b), links)
            (i === nothing || j === nothing) &&
                error("FreeSpace $name is not wired on both sides")
            src, _, τ1, φ1 = links[i]
            _, dst, τ2, φ2 = links[j]
            deleteat!(links, sort([i, j]))
            push!(links, (src, dst, τ1 + τ2 + r.τ, φ1 + φ2 + r.φ))
        end
    end

    src_of = Dict{Port, Int}()
    dst_of = Dict{Port, Int}()
    for (k, (src, dst, _, _)) in enumerate(links)
        haskey(src_of, src) && error("port $src sources two connections")
        haskey(dst_of, dst) && error("port $dst sinks two connections")
        src_of[src] = k
        dst_of[dst] = k
    end
    function edge(table, port)
        k = get(table, port, 0)
        k == 0 && error("optical port $(port[1]).$(port[2]) is not wired")
        return k
    end

    ext = [nameof(c) for c in comps if role[nameof(c)] isa Terminal]
    chan = Dict(name => e for (e, name) in enumerate(ext))
    K, E = length(links), length(ext)
    M1 = zeros(ComplexF64, K, K)
    M2 = zeros(ComplexF64, K, E)
    M3 = zeros(ComplexF64, E, K)
    M4 = zeros(ComplexF64, E, E)
    rot = [cis(l[4]) for l in links]    # folded propagator rotation per edge

    for (name, r) in role
        if r isa Terminal
            M2[edge(src_of, (name, :p1)), chan[name]] = 1
            k = edge(dst_of, (name, :p1))
            M3[chan[name], k] = rot[k]
        elseif r isa Scatterer
            for (o, i, c) in r.S
                k = edge(dst_of, (name, i))
                M1[edge(src_of, (name, o)), k] += c * rot[k]
            end
        end
    end

    x = [Symbol(s[1], '_', s[2], '→', d[1], '_', d[2]) for (s, d, _, _) in links]
    return (; M1, M2, M3, M4, delays = [l[3] for l in links], x, u = ext, y = copy(ext))
end
