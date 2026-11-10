# Lowering a composed optical model to the canonical transcendental scattering
# form: "rationalizing latest".
#
# For frequency-domain questions the model need not rationalize at all: every
# FreeSpace stays an exact delay `e^{-sτ}` on the edge it collapses into. The
# result is the data of Breadboard's `ScatteringSystem`:
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
field for field the contents of a Breadboard `ScatteringSystem`, so that

    G(s) = M3 * ((Diagonal(exp.(s .* delays)) - M1) \\ M2) + M4

is the exact transfer matrix from the external inputs `u` to the external
outputs `y`, FSR physics included. With Breadboard loaded,
`Breadboard.ScatteringSystem(to_scattering(model)...)` gives its solvers and
`pade` rationalization.

The lowering works on the scalar complex amplitude and reads only the
[`OpticalComponent`](@ref) metadata of each component: `x` indexes directed
edges, each [`Propagator`](@ref) (`FreeSpace`) collapses into an exact delay on
its edges, and each
[`Terminal`](@ref) (`Laser`, `Vacuum`, `Photodetector`, …) contributes one
input and one output channel named after it. Any `PropagationRep` is ignored;
the delays here are exact. Components scatter at their constructor values (a
`MovableMirror` at `x = 0`).

Subassemblies are flattened: a subsystem without `OpticalComponent` metadata is
searched for optical components, and its own ports may pass fields through to
them (`connect(p1.i, m1.p1.i)`, `connect(m1.p1.o, p1.o)`). Components are
named by their full path, e.g. `Symbol("cavity₊m1")`. Non-optical subsystems
and equations are ignored. Every optical port of every component must be
wired.
"""
function to_scattering(model::System)
    # Optical components (by full path) and the connectors each connect
    # equation joins, at every level of the hierarchy.
    role = Dict{Symbol, Any}()
    order = Symbol[]
    joined = Vector{Vector{String}}()
    function walk(sys, prefix)
        for eq in ModelingToolkit.get_eqs(sys)
            rhs = SymbolicUtils.unwrap_const(Symbolics.unwrap(eq.rhs))
            rhs isa ModelingToolkit.Connection || continue
            push!(joined, [prefix * string(nameof(c)) for c in rhs.systems])
        end
        for sub in ModelingToolkit.get_systems(sys)
            path = prefix * string(nameof(sub))
            r = SymbolicUtils.getmetadata(sub, OpticalComponent, nothing)
            if r === nothing
                walk(sub, path * "₊")
            else
                role[Symbol(path)] = r
                push!(order, Symbol(path))
            end
        end
        return
    end
    walk(model, "")

    # Connected connectors carry one field. Union them; each class must hold
    # exactly one component `o` (where the field is made) and one component
    # `i` (where it arrives), plus any subassembly ports it passes through.
    parent = Dict{String, String}()
    root(c) = (p = get!(parent, c, c); p == c ? c : (parent[c] = root(p)))
    for cs in joined, c in cs[2:end]
        parent[root(c)] = root(cs[1])
    end
    Port = Tuple{Symbol, Symbol}
    ends = Dict{String, Dict{String, Vector{Port}}}()
    for c in keys(parent)
        parts = split(c, '₊')
        length(parts) ≥ 3 || continue
        owner = Symbol(join(parts[1:(end - 2)], '₊'))
        haskey(role, owner) && parts[end] in ("i", "o") || continue
        class = get!(() -> Dict("i" => Port[], "o" => Port[]), ends, root(c))
        push!(class[parts[end]], (owner, Symbol(parts[end - 1])))
    end

    # Directed links (src, dst, delay) between component ports.
    links = Tuple{Port, Port, Float64}[]
    for class in values(ends)
        o, i = class["o"], class["i"]
        length(o) == 1 && length(i) == 1 || error(
            "a field must leave one port and arrive at one port; found it leaving " *
                "$(join(o, ", ", " and ")) and arriving at $(join(i, ", ", " and "))"
        )
        push!(links, (only(o), only(i), 0.0))
    end
    sort!(links)

    # Collapse propagators: the link into one side and the link out of the
    # other fuse into a single edge carrying the delay.
    for name in order
        r = role[name]
        r isa Propagator || continue
        for (a, b) in ((:p1, :p2), (:p2, :p1))
            i = findfirst(l -> l[2] == (name, a), links)
            j = findfirst(l -> l[1] == (name, b), links)
            (i === nothing || j === nothing) &&
                error("FreeSpace $name is not wired on both sides")
            src, _, τ1 = links[i]
            _, dst, τ2 = links[j]
            deleteat!(links, sort([i, j]))
            push!(links, (src, dst, τ1 + τ2 + r.τ))
        end
    end

    src_of = Dict(l[1] => k for (k, l) in enumerate(links))
    dst_of = Dict(l[2] => k for (k, l) in enumerate(links))
    function edge(table, port)
        k = get(table, port, 0)
        k == 0 && error("optical port $(port[1]).$(port[2]) is not wired")
        return k
    end

    ext = [name for name in order if role[name] isa Terminal]
    chan = Dict(name => e for (e, name) in enumerate(ext))
    K, E = length(links), length(ext)
    M1 = zeros(ComplexF64, K, K)
    M2 = zeros(ComplexF64, K, E)
    M3 = zeros(ComplexF64, E, K)
    M4 = zeros(ComplexF64, E, E)

    for name in order
        r = role[name]
        if r isa Terminal
            M2[edge(src_of, (name, :p1)), chan[name]] = 1
            M3[chan[name], edge(dst_of, (name, :p1))] = 1
        elseif r isa Scatterer
            for (o, i, c) in r.S
                M1[edge(src_of, (name, o)), edge(dst_of, (name, i))] += c
            end
        end
    end

    x = [Symbol(s[1], '_', s[2], '→', d[1], '_', d[2]) for (s, d, _) in links]
    return (; M1, M2, M3, M4, delays = [l[3] for l in links], x, u = ext, y = copy(ext))
end
