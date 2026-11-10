# KatScript: one model file, read by both OpticalModelingToolkit and Finesse.jl.
#
# A `.kat` file is the single source of truth for a model. Finesse.jl parses the
# whole language and solves the frequency response. OMT reads the subset below
# and builds the same network as a time-domain DAE. The two readings are checked
# against each other in `test/finesse/`.
#
# Supported elements (a subset of KatScript):
#
#     laser  name P=<W>
#     mirror name R=<power> [T=<power>] [L=<power>] [phi=<deg>]
#     space  name a.pN b.pN L=<m>
#     pendulum name m.mech mass=<kg> fz=<Hz> [Qz=<Q>]
#
# Analysis lines (`fsig`, `frequency_response`, `xaxis`, `noxaxis`, `readout_dc`)
# are accepted and ignored: they say what Finesse should compute, not what the
# network is. Anything else is an error.
#
# Conventions, as established by the Finesse cross-check:
#
#   * A Finesse mirror is a turned-around OMT mirror: Finesse `pN` is OMT
#     `p(3-N)`. Unconnected optical ports are terminated with `Vacuum`, which is
#     what Finesse does with an open node.
#   * Finesse tuning `phi` [deg] is OMT tuning `ϕ` [rad] with the opposite sign.
#   * A `space` is a pure delay. Its carrier phase is whole wavelengths, so only
#     the sidebands see `L`, as in `FreeSpace`.
#   * A pendulum is `M z̈ + (M ω₀/Q) ż + M ω₀² z = F_rad + F_z`, with `ω₀ = 2π fz`.
#     `z` is the displacement of the mirror it suspends and `F_z` the external
#     force on it. It is the same transfer function as Finesse's `pendulum`.
#     The mirror's tuning follows `-z`, which is the sign Finesse uses.

const KAT_ELEMENTS = (:laser, :mirror, :space, :pendulum)
const KAT_ANALYSIS = (:fsig, :frequency_response, :xaxis, :noxaxis, :readout_dc)
const KAT_WAVELENGTH = 1064.0e-9

"""
    KatElement(kind, name, nodes, params)

One element line of a KatScript model. `nodes` are the port references on the
line, as `(element, port)` pairs; `params` are its numeric `key=value` arguments.
"""
struct KatElement
    kind::Symbol
    name::Symbol
    nodes::Vector{Tuple{Symbol, Symbol}}
    params::Dict{Symbol, Float64}
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
        ## analysis calls such as `fsig(1)` or `frequency_response(...)`
        occursin('(', line) && continue
        toks = split(line)
        kind = Symbol(toks[1])
        kind in KAT_ANALYSIS && continue
        kind in KAT_ELEMENTS ||
            throw(ArgumentError("KatScript: unsupported element `$kind` in `$line`"))
        length(toks) >= 2 || throw(ArgumentError("KatScript: no name in `$line`"))
        nodes = Tuple{Symbol, Symbol}[]
        params = Dict{Symbol, Float64}()
        for tok in toks[3:end]
            if occursin('=', tok)
                key, value = split(tok, '='; limit = 2)
                v = tryparse(Float64, value)
                v === nothing &&
                    throw(ArgumentError("KatScript: `$key` must be a plain number in `$line`"))
                params[Symbol(key)] = v
            elseif occursin('.', tok)
                elem, port = split(tok, '.'; limit = 2)
                push!(nodes, (Symbol(elem), Symbol(port)))
            else
                throw(ArgumentError("KatScript: cannot read `$tok` in `$line`"))
            end
        end
        push!(elements, KatElement(kind, Symbol(toks[2]), nodes, params))
    end
    return elements
end

"""
    build_kat(text; drive = nothing, rep = Algebraic(), x_rest = 0.0)
        -> (plant, handles)

Build the time-domain model of a KatScript network. The returned `plant` is a
`System` with the pendulum's force `F_z` as an input (for `linearize_proper` and
the frequency response) and its displacement `z` as the output.

With `drive = (A, f)` the force is instead `F_z = A sin(2π f t)` and `F_z` is not an
input, which is what a time-domain simulation needs.

`x_rest` is the rest position of the suspension. Set it so the static radiation
force is balanced at the operating point, as `test/test_optomechanics.jl` does.

`handles` is a `NamedTuple` of the variables and components a caller needs:
`z`, `v`, `F_z`, and `comps`, a `Dict` of the component systems by name.
"""
function build_kat(text::AbstractString; drive = nothing, rep::PropagationRep = Algebraic(), x_rest = 0.0)
    elements = parse_kat(text)
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
    for e in elements
        if e.kind === :laser
            haskey(e.params, :P) || throw(ArgumentError("KatScript: laser `$(e.name)` needs P="))
            comps[e.name] = Laser(; name = e.name, A = sqrt(e.params[:P]))
        elseif e.kind === :mirror
            haskey(e.params, :R) || throw(ArgumentError("KatScript: mirror `$(e.name)` needs R="))
            R = e.params[:R]
            L = get(e.params, :L, 0.0)
            T = get(e.params, :T, 1 - R - L)
            ϕ = -deg2rad(get(e.params, :phi, 0.0))
            comps[e.name] = haskey(suspended, e.name) ?
                MovableMirror(; name = e.name, R, T, ϕ, λ = KAT_WAVELENGTH) :
                Mirror(; name = e.name, R, T, ϕ)
        elseif e.kind === :space
            haskey(e.params, :L) || throw(ArgumentError("KatScript: space `$(e.name)` needs L="))
            comps[e.name] = FreeSpace(; name = e.name, L = e.params[:L], rep)
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
        ports = e.kind === :laser ? (:p1,) : (:p1, :p2)
        for port in ports
            (e.name, port) in connected && continue
            vac = Vacuum(; name = Symbol(e.name, :_, port, :_vac))
            push!(systems, vac)
            append!(eqs, link(vac.p1, omt_port(e.name, port)))
        end
    end

    if isempty(pendulums)
        @named plant = System(eqs, t; systems)
        return plant, (; comps)
    end

    p = only(pendulums)
    m = p.nodes[1][1]
    mirror = comps[m]
    M = p.params[:mass]
    Q = get(p.params, :Qz, 1000.0)
    ω0 = 2π * p.params[:fz]
    k = M * ω0^2
    if drive === nothing
        @variables z(t) = 0.0 v(t) = 0.0 F_z(t) = 0.0 [input = true]
    else
        A, f = drive
        @variables z(t) = 0.0 v(t) = 0.0 F_z(t) = 0.0
        push!(eqs, F_z ~ A * sin(2π * f * t))
    end
    append!(
        eqs, [
            ## Finesse's displacement `z` moves the tuning with the opposite sign to
            ## OMT's `x`. Found by the suspended-cavity cross-check at ±0.5 deg, where the
            ## optical spring is nonzero; it has no effect at resonance.
            mirror.x ~ -z
            D(z) ~ v
            M * D(v) ~ -k * (z - x_rest) - (M * ω0 / Q) * v + mirror.F + F_z
        ]
    )
    @named plant = System(eqs, t; systems)
    return plant, (; z, v, F_z, comps)
end
