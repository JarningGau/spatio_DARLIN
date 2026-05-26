#!/usr/bin/env bash
set -euo pipefail

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

command -v snakemake >/dev/null || fail "snakemake not found (activate the spatio_darlin conda env)"
command -v cutadapt >/dev/null || fail "cutadapt not found"
command -v umi_tools >/dev/null || fail "umi_tools not found"
python -c "import darlin_core, spatio_darlin" || fail "darlin_core or spatio_darlin not importable"
command -v BSTMatrix >/dev/null || fail "BSTMatrix not on PATH (see README step 1)"

echo "Analysis environment is ready!"
