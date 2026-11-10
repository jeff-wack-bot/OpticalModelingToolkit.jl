# The shared models in `models/` and the Finesse.jl responses kept as references.
#
# `models/reference/<model>.txt` holds Finesse's frequency responses of
# `models/<model>.kat`, so the main test suite and the docs can compare OMT
# against Finesse without Finesse installed. `reference.jl` writes them;
# `test_kat_models.jl` checks they are current.
#
# `locked_cavity.kat` wires a servo into a mirror's force node, which needs
# Finesse.jl's `fix/wire-node-ends` (submitted upstream); with an older
# Finesse.jl, its reference and cross-check fail at compile time.

const MODELS_DIR = joinpath(@__DIR__, "..", "..", "models")

"Frequencies [Hz] of every stored reference: 31 points, 1 Hz to 1 kHz."
const REFERENCE_FREQS = exp10.(range(0, 3; length = 31))

"Each model, and the Finesse `(input, output)` node pairs stored for it."
const REFERENCES = [
    "suspended_fp" => [("m2.mech.F_z", "m2.mech.z")],
    "locked_cavity" => [("m2.mech.F_z", "m2.mech.z"), ("m2.mech.F_z", "trans.DC.o")],
]

kat_source(model) = read(joinpath(MODELS_DIR, "$model.kat"), String)
reference_path(model) = joinpath(MODELS_DIR, "reference", "$model.txt")

"Finesse's frequency response from node `input` to node `output` of `src` at `fs` [Hz]."
function finesse_response(src, fs, input, output)
    freqs = "[" * join(fs, ", ") * "]"
    sol = finesse_run(src * "\nfsig(1)\nfrequency_response($freqs, $input, $output)\n")
    return vec(sol.data)
end

"Every stored response of `model`, freshly computed by Finesse: `pairs => columns`."
function finesse_references(model, pairs)
    src = kat_source(model)
    return [finesse_response(src, REFERENCE_FREQS, a, b) for (a, b) in pairs]
end

function write_reference(model, pairs, responses)
    open(reference_path(model), "w") do io
        println(io, "# Finesse.jl frequency responses of models/$model.kat.")
        println(io, "# Written by test/finesse/reference.jl; test/finesse checks it is current.")
        println(io, "# responses: ", join(("$a->$b" for (a, b) in pairs), " "))
        println(io, "# columns: f [Hz], then the real and imaginary part of each response")
        for (i, f) in enumerate(REFERENCE_FREQS)
            cols = [repr(f); [repr(x) for r in responses for x in (real(r[i]), imag(r[i]))]]
            println(io, join(cols, " "))
        end
    end
    return nothing
end

"Read a stored reference: `(freqs, responses)`, each response a complex vector."
function read_reference(model)
    rows = [parse.(Float64, split(l)) for l in eachline(reference_path(model)) if !startswith(l, '#')]
    fs = [r[1] for r in rows]
    n = (length(rows[1]) - 1) ÷ 2
    return fs, [[complex(r[2k], r[2k + 1]) for r in rows] for k in 1:n]
end
