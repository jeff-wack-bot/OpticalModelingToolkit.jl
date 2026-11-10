# KatScript: one model file, read by both OpticalModelingToolkit and Finesse.jl.
#
# A `.kat` file is the single source of truth for a model. Finesse.jl parses the
# whole language and solves the frequency response. OMT reads the subset below
# and builds the same network as a time-domain DAE, so the nonlinear and
# transient questions Finesse cannot ask (lock acquisition, large excursions)
# are asked of the same model. The two readings are checked against each other
# in `test/finesse/`.
#
# Supported lines (a subset of KatScript; numbers must be plain literals):
#
#     laser      name P=<W>
#     mirror     name R=<power> [T=<power>] [L=<power>] [phi=<deg>]
#     space      name a.pN b.pN L=<m>
#     pendulum   name m.mech mass=<kg> fz=<Hz> [Qz=<Q>]
#     readout_dc name m.pN.(i|o)
#     zpk        name [zeros] [poles] [k=<k>] [gain=<g>]     (roots in rad/s)
#     wire       name from to [gain=<g>]
#     link(from, to, ...)                                      (signal nodes only)
#
# Wires carry signals: from a readout `rd.DC.o` or a filter output `f.p2.o`, to a
# filter input `f.p1.i` or a suspended mirror's force `m.mech.F_z`. Wires into
# the same node add.
#
# Analysis lines (`fsig`, `frequency_response`, `xaxis`, `noxaxis`, and any line
# with a call) are accepted and ignored: they say what Finesse should compute,
# not what the network is. Anything else is an error.
#
# Conventions, as established by the Finesse cross-check:
#
#   * A Finesse mirror is a turned-around OMT mirror: Finesse `pN` is OMT
#     `p(3-N)`. Unconnected optical ports are terminated with `Vacuum`, which is
#     what Finesse does with an open node.
#   * Finesse tuning `phi` [deg] is OMT tuning `ϕ` [rad] with the opposite sign,
#     and a suspended mirror's tuning follows `-z`.
#   * A `space` is a pure delay. Its carrier phase is whole wavelengths, so only
#     the sidebands see `L`, as in `FreeSpace`.
#   * A pendulum is `M z̈ + (M ω₀/Q) ż + M ω₀² (z - z_rest) = F`, with `ω₀ = 2π fz`
#     and `F` the radiation force plus every force wired or applied to it: the
#     transfer function of Finesse's `pendulum`.
#   * Finesse linearizes about the DC state with every mirror at its tuning. So
#     OMT puts the model in that state: `z_rest` balances the static radiation
#     force at `z = 0`, and each readout's DC power is subtracted from it, so
#     the error signals and filter outputs are zero there.

const KAT_ELEMENTS = (:laser, :mirror, :space, :pendulum, :readout_dc, :zpk, :wire)
const KAT_ANALYSIS = (:fsig, :frequency_response, :xaxis, :noxaxis)
const KAT_WAVELENGTH = 1064.0e-9

"""
    KatElement(kind, name, nodes, params, lists)

One element line of a KatScript model. `nodes` are the port references on the
line, as `(element, rest)` pairs with `rest` the remainder of the dotted path
(`:p1`, `Symbol("p1.o")`, `Symbol("mech.F_z")`, ...); `params` are its numeric
`key=value` arguments and `lists` its bracketed number lists, positional ones
under `:_1`, `:_2`, ... and named ones under their key.
"""
struct KatElement
    kind::Symbol
    name::Symbol
    nodes::Vector{Tuple{Symbol, Symbol}}
    params::Dict{Symbol, Float64}
    lists::Dict{Symbol, Vector{Float64}}
end

function parse_kat_number(s, line)
    v = tryparse(Float64, s)
    v === nothing && throw(ArgumentError("KatScript: `$s` must be a plain number in `$line`"))
    return v
end

function parse_kat_list(s, line)
    body = strip(s[2:(end - 1)])
    isempty(body) && return Float64[]
    return [parse_kat_number(strip(x), line) for x in split(body, ',')]
end

"""
    parse_kat(text) -> Vector{KatElement}

Parse the supported subset of KatScript. Comments (`#`), blank lines and analysis
lines are skipped. Unknown elements, malformed arguments and non-numeric values
are errors.
"""
function parse_kat(text::AbstractString)
    elements = KatElement[]
    for raw in eachline(IOBuffer(text))
        line = strip(first(split(raw, '#'; limit = 2)))
        isempty(line) && continue
        ## `link(a, b, c)` between signal nodes is wires a→b→c; it must not be skipped
        ## with the analysis calls, or a servo would silently lose a connection
        if startswith(line, "link(")
            append!(elements, kat_link(line, length(elements)))
            continue
        end
        ## analysis calls such as `fsig(1)` or `frequency_response(...)`
        occursin('(', line) && continue
        ## a list may contain spaces; close them up so it is one token
        line = replace(line, r"\[[^\]]*\]" => m -> replace(m, r"\s+" => ""))
        toks = split(line)
        kind = Symbol(toks[1])
        kind in KAT_ANALYSIS && continue
        kind in KAT_ELEMENTS ||
            throw(ArgumentError("KatScript: unsupported element `$kind` in `$line`"))
        length(toks) >= 2 || throw(ArgumentError("KatScript: no name in `$line`"))
        nodes = Tuple{Symbol, Symbol}[]
        params = Dict{Symbol, Float64}()
        lists = Dict{Symbol, Vector{Float64}}()
        npos = 0
        for tok in toks[3:end]
            if startswith(tok, '[')
                lists[Symbol(:_, npos += 1)] = parse_kat_list(tok, line)
            elseif occursin('=', tok)
                key, value = split(tok, '='; limit = 2)
                if startswith(value, '[')
                    lists[Symbol(key)] = parse_kat_list(value, line)
                else
                    params[Symbol(key)] = parse_kat_number(value, line)
                end
            elseif occursin('.', tok)
                elem, rest = split(tok, '.'; limit = 2)
                push!(nodes, (Symbol(elem), Symbol(rest)))
            else
                throw(ArgumentError("KatScript: cannot read `$tok` in `$line`"))
            end
        end
        push!(elements, KatElement(kind, Symbol(toks[2]), nodes, params, lists))
    end
    return elements
end

## `link(a, b, ...)` with every argument a dotted signal node, as wires between
## consecutive arguments. An optical link (between elements or ports) is refused:
## OMT models a connection between optical ports as a `space`.
function kat_link(line, n)
    m = match(r"^link\((.*)\)$", line)
    m === nothing && throw(ArgumentError("KatScript: cannot read `$line`"))
    args = [strip(a) for a in split(m[1], ',') if !occursin('=', a)]
    nodes = map(args) do a
        parts = split(a, '.')
        length(parts) >= 3 || throw(
            ArgumentError(
                "KatScript: `link` is supported only between signal nodes (like `rd.DC.o`); " *
                    "connect optical ports with a `space`, in `$line`"
            )
        )
        (Symbol(parts[1]), Symbol(join(parts[2:end], '.')))
    end
    length(nodes) >= 2 || throw(ArgumentError("KatScript: `link` needs two nodes in `$line`"))
    return [
        KatElement(:wire, Symbol(:link, n + i), [nodes[i], nodes[i + 1]], Dict{Symbol, Float64}(), Dict{Symbol, Vector{Float64}}())
            for i in 1:(length(nodes) - 1)
    ]
end

kat_require(e, key) = haskey(e.params, key) ? e.params[key] :
    throw(ArgumentError("KatScript: $(e.kind) `$(e.name)` needs $key="))

"""
    build_kat(text; drive = nothing, rep = Algebraic(), operating_point = true)
        -> (plant, handles)

Build the time-domain model of a KatScript network.

The pendulum's force `F_z` is an input and its displacement `z` an output, ready
for [`linearize_proper`](@ref). With `drive = (A, f)` the force is instead
`F_z = A sin(2π f t)` and is not an input, which is what a time-domain simulation
needs. `A` and `f` are the parameters `handles.drive.drive_A` and `drive_f`, so
one compiled system serves for every drive; see [`steady_state`](@ref).

With `operating_point = true` the model is put in the state Finesse linearizes
about (see the conventions at the top of `src/katscript.jl`), found by solving the
static network once. `rep` chooses how every `space` represents its delay; the
optical spring's damping comes from the cavity's storage time, so a detuned
suspended cavity needs `PadeDelay`, not the adiabatic `Algebraic`.

`handles` holds what a caller needs: `z`, `v`, `F_z`, the component systems
`comps`, the readout powers `readouts` and the filters `filters` (each a `Dict` by
name), and the operating point `op` (`z_rest` and the readout `offsets`).
"""
function build_kat(text::AbstractString; drive = nothing, rep::PropagationRep = Algebraic(), operating_point = true)
    elements = parse_kat(text)
    op = (; z_rest = 0.0, offsets = Dict{Symbol, Float64}())
    if operating_point
        op = kat_operating_point(elements)
    end
    return kat_assemble(elements; drive, rep, op, static = false)
end

"""
    kat_operating_point(elements) -> (; z_rest, offsets)

Solve the static network, with any suspended mirror held at its tuning (`z = 0`),
and return the suspension rest position that balances the radiation force there
and the DC power of each readout.
"""
function kat_operating_point(elements::Vector{KatElement})
    plant, h = kat_assemble(elements; drive = nothing, rep = Algebraic(), op = nothing, static = true)
    sys = mtkcompile(plant)
    prob = NonlinearProblem(sys, Dict(u => 0.0 for u in unknowns(sys)))
    sol = SciMLBase.solve(prob, SimpleNewtonRaphson(); abstol = 1.0e-14, reltol = 1.0e-14)
    SciMLBase.successful_retcode(sol) ||
        error("KatScript: the static network did not solve ($(sol.retcode))")
    offsets = Dict{Symbol, Float64}(name => sol[P] for (name, P) in h.readouts)
    z_rest = 0.0
    if h.suspension !== nothing
        F0 = sol[h.suspension.mirror.F]
        z_rest = -F0 / h.suspension.k
    end
    return (; z_rest, offsets)
end

## The model. `static = true` holds a suspended mirror at its tuning and leaves out
## the suspension and the electronics, which is the network `kat_operating_point`
## solves; `op` is then unused.
function kat_assemble(elements; drive, rep, op, static)
    pendulums = filter(e -> e.kind === :pendulum, elements)
    length(pendulums) <= 1 ||
        throw(ArgumentError("build_kat supports one pendulum per model; found $(length(pendulums))"))
    suspended = Dict{Symbol, KatElement}()
    for p in pendulums
        length(p.nodes) == 1 && p.nodes[1][2] === :mech ||
            throw(ArgumentError("KatScript: pendulum `$(p.name)` must attach to `<mirror>.mech`"))
        suspended[p.nodes[1][1]] = p
    end

    comps = Dict{Symbol, Any}()
    filters = Dict{Symbol, Any}()
    for e in elements
        if e.kind === :laser
            comps[e.name] = Laser(; name = e.name, A = sqrt(kat_require(e, :P)))
        elseif e.kind === :mirror
            R = kat_require(e, :R)
            T = get(e.params, :T, 1 - R - get(e.params, :L, 0.0))
            ϕ = -deg2rad(get(e.params, :phi, 0.0))
            comps[e.name] = haskey(suspended, e.name) ?
                MovableMirror(; name = e.name, R, T, ϕ, λ = KAT_WAVELENGTH) :
                Mirror(; name = e.name, R, T, ϕ)
        elseif e.kind === :space
            comps[e.name] = FreeSpace(; name = e.name, L = kat_require(e, :L), rep)
        elseif e.kind === :zpk && !static
            z = get(e.lists, :z, get(e.lists, :_1, Float64[]))
            p = get(e.lists, :p, get(e.lists, :_2, Float64[]))
            k = haskey(e.params, :k) ? e.params[:k] : nothing
            filters[e.name] = ZPK(; name = e.name, z, p, k, gain = get(e.params, :gain, 1.0))
        end
    end

    kinds = Dict(e.name => e.kind for e in elements)

    ## the OMT port that Finesse port `elem.port` refers to. A mirror is turned
    ## around (Finesse p1 is OMT p2); a laser has one port, p1.
    function omt_port(elem::Symbol, port::Symbol)
        haskey(kinds, elem) || throw(ArgumentError("KatScript: unknown element `$elem`"))
        kind = kinds[elem]
        if kind === :laser && port === :p1
            return comps[elem].p1
        elseif kind === :mirror && port === :p1
            return comps[elem].p2
        elseif kind === :mirror && port === :p2
            return comps[elem].p1
        end
        throw(ArgumentError("KatScript: `$elem.$port` is not an optical port"))
    end

    eqs = Equation[]
    connected = Set{Tuple{Symbol, Symbol}}()
    systems = Any[values(comps)...]
    for e in elements
        e.kind === :space || continue
        length(e.nodes) == 2 ||
            throw(ArgumentError("KatScript: space `$(e.name)` needs two port references"))
        (a, b) = e.nodes
        append!(eqs, link(comps[e.name].p1, omt_port(a...)))
        append!(eqs, link(comps[e.name].p2, omt_port(b...)))
        push!(connected, a, b)
    end

    ## terminate every optical port that no line connects, as Finesse leaves it dark
    for e in elements
        e.kind in (:mirror, :laser) || continue
        for port in (e.kind === :laser ? (:p1,) : (:p1, :p2))
            (e.name, port) in connected && continue
            vac = Vacuum(; name = Symbol(e.name, :_, port, :_vac))
            push!(systems, vac)
            append!(eqs, link(vac.p1, omt_port(e.name, port)))
        end
    end

    ## readouts: the power at a node, less its DC value at the operating point
    readouts = Dict{Symbol, Any}()
    for e in elements
        e.kind === :readout_dc || continue
        length(e.nodes) == 1 ||
            throw(ArgumentError("KatScript: readout_dc `$(e.name)` needs one optical node"))
        elem, rest = e.nodes[1]
        parts = split(String(rest), '.')
        length(parts) == 2 && parts[2] in ("i", "o") ||
            throw(ArgumentError("KatScript: readout_dc `$(e.name)` needs a node like `m.p1.o`"))
        port = omt_port(elem, Symbol(parts[1]))
        conn = parts[2] == "o" ? port.o : port.i
        P = only(@variables $(Symbol(e.name, :_DC))(t))
        offset = static ? 0.0 : get(op.offsets, e.name, 0.0)
        push!(eqs, P ~ abs2(field(conn)) - offset)
        readouts[e.name] = P
    end

    ## signal nodes that wires read from and add into
    inflow = Dict{Tuple{Symbol, Symbol}, Any}()          # (element, node) => sum
    function source(elem, rest)
        rest === Symbol("DC.o") && haskey(readouts, elem) && return readouts[elem]
        rest === Symbol("p2.o") && haskey(filters, elem) && return filters[elem].y
        throw(ArgumentError("KatScript: `$elem.$rest` is not a signal output"))
    end
    function sink!(elem, rest, value)
        ok = (rest === Symbol("p1.i") && haskey(filters, elem)) ||
            (rest === Symbol("mech.F_z") && haskey(suspended, elem))
        ok || throw(ArgumentError("KatScript: `$elem.$rest` is not a signal input"))
        key = (elem, rest)
        inflow[key] = get(inflow, key, 0) + value
        return nothing
    end
    if !static
        for e in elements
            e.kind === :wire || continue
            length(e.nodes) == 2 ||
                throw(ArgumentError("KatScript: wire `$(e.name)` needs a source and a destination"))
            get(e.params, :delay, 0.0) == 0 ||
                throw(ArgumentError("KatScript: wire `$(e.name)`: delays are not supported"))
            (from, to) = e.nodes
            sink!(to..., get(e.params, :gain, 1.0) * source(from...))
        end
        for (name, f) in filters
            push!(systems, f)
            push!(eqs, f.u ~ get(inflow, (name, Symbol("p1.i")), 0.0))
        end
    end

    handles = (; comps, readouts, filters)
    if isempty(pendulums)
        @named plant = System(eqs, t; systems)
        return plant, (; handles..., z = nothing, v = nothing, F_z = nothing, suspension = nothing, op, drive = ())
    end

    p = only(pendulums)
    m = p.nodes[1][1]
    mirror = comps[m]
    M = kat_require(p, :mass)
    Q = get(p.params, :Qz, 1000.0)
    ω0 = 2π * kat_require(p, :fz)
    k = M * ω0^2
    suspension = (; mirror, k)
    if static
        push!(eqs, mirror.x ~ 0)
        @named plant = System(eqs, t; systems)
        return plant, (; handles..., z = nothing, v = nothing, F_z = nothing, suspension, op, drive = ())
    end

    if drive === nothing
        @variables z(t) = 0.0 v(t) = 0.0 F_z(t) = 0.0 [input = true]
        drive_pars = ()
    else
        A, f = drive
        @variables z(t) = 0.0 v(t) = 0.0 F_z(t) = 0.0
        @parameters drive_A = A [description = "Drive amplitude [N]"] drive_f = f [description = "Drive frequency [Hz]"]
        push!(eqs, F_z ~ drive_A * sin(2π * drive_f * t))
        drive_pars = (; drive_A, drive_f)
    end
    F_wired = get(inflow, (m, Symbol("mech.F_z")), 0.0)
    append!(
        eqs, [
            mirror.x ~ -z
            D(z) ~ v
            M * D(v) ~ -k * (z - op.z_rest) - (M * ω0 / Q) * v + mirror.F + F_z + F_wired
        ]
    )
    @named plant = System(eqs, t; systems)
    return plant, (; handles..., z, v, F_z, suspension, op, drive = drive_pars)
end
