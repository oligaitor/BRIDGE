#!/bin/bash
# Create BRIDGE's own Python environment (a venv with the Python packages).
#
#   ./setup_env.sh                 # environment in .venv next to this file
#   ./setup_env.sh --slim          # ...as few files as possible (see below)
#   VENV=/somewhere/else ./setup_env.sh
#   ./setup_env.sh --pinned        # exact versions from requirements.txt
#   ./setup_env.sh --dev           # also pytest, for running the tests
#
# Needs Python 3.10+ as `python3` (or set PYTHON=/path/to/python3.1x). On a
# cluster this usually means loading a Python module first.
#
# File count: many HPC project quotas limit the number of files as well as the
# space, and a normal install is about 8,000 small files. --slim cuts that to
# about 3,200 by not writing .pyc bytecode and deleting the packages' own test
# suites, neither of which BRIDGE needs.
#
# Afterwards, every session:
#   source <venv>/bin/activate && source tools.env
#   bridge check
#
# The venv holds only Python and the Python packages. The binaries the arms
# call (kallisto, mash, minimap2, bwa-mem2, samtools, BLAST, viral_consensus,
# CliqueSNV + Java, KMA) come from tools.env, or from a conda environment
# built from environment.yml.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
py=${PYTHON:-python3}
venv=${VENV:-$here/.venv}
pinned=0
dev=0
slim=0
for arg in "$@"; do
    case "$arg" in
        --pinned) pinned=1 ;;
        --dev) dev=1 ;;
        --slim) slim=1 ;;
        --venv=*) venv=${arg#--venv=} ;;
        *) echo "unknown option: $arg" >&2; exit 2 ;;
    esac
done
if (( slim && dev )); then
    echo "ERROR: --slim removes the test suites; --dev installs them. Pick one." >&2
    exit 2
fi

if ! "$py" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' 2>/dev/null; then
    found=$("$py" -c 'import sys; print(".".join(map(str, sys.version_info[:3])))' 2>/dev/null || echo "none")
    echo "ERROR: need Python 3.10+ on PATH as '$py' (found: $found)." >&2
    echo "       On a cluster, load a Python 3.10+ module first, then re-run." >&2
    echo "       Or set PYTHON=/path/to/python3.1x" >&2
    exit 1
fi

quota_help() {
    echo "       'Disk quota exceeded' means that filesystem is out of space OR out of" >&2
    echo "       FILES (check your site's quota tool for a file count)." >&2
    echo "       Then either build a slim environment:" >&2
    echo "         rm -rf \"$venv\" && $0 --slim" >&2
    echo "       or put it somewhere with more room:" >&2
    echo "         rm -rf \"$venv\" && VENV=/path/with/more/room/bridge_venv $0" >&2
}

if [[ ! -x "$venv/bin/python" ]]; then
    echo "Creating $venv with $("$py" -c 'import sys; print(sys.version.split()[0])') ..."
    if ! "$py" -m venv "$venv"; then
        echo "ERROR: could not create $venv." >&2
        quota_help
        exit 1
    fi
fi

install() {
    local options=(--quiet --no-cache-dir)
    (( slim )) && options+=(--no-compile)
    if ! "$venv/bin/python" -m pip install "${options[@]}" "$@"; then
        echo "ERROR: pip install failed (see above)." >&2
        quota_help
        exit 1
    fi
}

install --upgrade pip
if (( pinned )); then
    install -r "$here/requirements.txt"
fi
if (( dev )); then
    install -e "$here[dev]"
else
    install -e "$here"
fi

if (( slim )); then
    # The packages' own test suites: tens of megabytes and thousands of files
    # that nothing here imports. Only directories named exactly "tests" go -
    # numpy.testing and pandas._testing are real modules and stay.
    echo "Trimming test suites and bytecode ..."
    find "$venv" -type d -name tests -prune -exec rm -rf {} + 2>/dev/null || true
    find "$venv" -type d -name __pycache__ -prune -exec rm -rf {} + 2>/dev/null || true
fi

echo
echo "Installed: $("$venv/bin/bridge" --version)"
echo "Environment: $venv"
echo "             $(find "$venv" -type f | wc -l) files, $(du -sh "$venv" 2>/dev/null | cut -f1)"
echo
echo "Each session, before running anything:"
echo "  source $venv/bin/activate && source $here/tools.env"
echo "  bridge check"
