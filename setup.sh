#!/usr/bin/env bash
#
# setup.sh — bootstrap a fresh checkout of OpticalModelingToolkit.jl for
# development.
#
# What this does:
#   1. Clones finesse3, ModelingToolkit.jl, and ModelingToolkitStandardLibrary.jl
#      into context/ for offline reading and grep-able reference. Idempotent:
#      already-present clones are left alone.
#   2. Reminds you that Python setup is automatic — CondaPkg.toml at the project
#      root tells PythonCall.jl how to materialize a project-local Python env
#      with finesse installed via pip. Nothing to install by hand.
#
# What this does NOT do:
#   - It does not run `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
#     for you. Run that yourself once setup.sh finishes; on first load of
#     PythonCall, CondaPkg will create .CondaPkg/ and download Python + finesse.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
context_dir="${repo_root}/context"
mkdir -p "${context_dir}"

clone_if_missing() {
    local name="$1"
    local url="$2"
    local target="${context_dir}/${name}"
    if [[ -d "${target}/.git" ]]; then
        echo "[setup] ${name} already present at ${target} — skipping"
    else
        echo "[setup] cloning ${name} from ${url}"
        git clone --depth 1 "${url}" "${target}"
    fi
}

clone_if_missing finesse3                          https://git.ligo.org/finesse/finesse3.git
clone_if_missing ModelingToolkit.jl                https://github.com/SciML/ModelingToolkit.jl.git
clone_if_missing ModelingToolkitStandardLibrary.jl https://github.com/SciML/ModelingToolkitStandardLibrary.jl.git
clone_if_missing sflu_optomechanics                git@git.ligo.org:jeffrey.wack/sflu-optomechanics.git
clone_if_missing pygwinc                           git@git.ligo.org:jeffrey.wack/pygwinc.git

cat <<'EOF'

[setup] context/ is ready.

Next step:

    julia --project=. -e 'using Pkg; Pkg.instantiate()'

The first time anything in this project does `using PythonCall`, CondaPkg will
read CondaPkg.toml and download Python + numpy + scipy + finesse into
.CondaPkg/. Expect that first load to take a few minutes. After that, no
further Python setup is required.

EOF
