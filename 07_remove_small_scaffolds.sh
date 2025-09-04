#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: bash 07_remove_small_scaffolds.sh <input.fq> <min_length>" >&2
  exit 1
fi

IN_FQ=$1
MINLEN=$2
OUT_FQ="${IN_FQ%.fq}.min${MINLEN}.fq"

# Require seqkit
if ! command -v seqkit &>/dev/null; then
  echo "Error: seqkit not found in PATH. Install with conda (bioconda::seqkit) or download binary." >&2
  exit 1
fi

# Filter sequences by length
seqkit seq -m "${MINLEN}" "${IN_FQ}" > "${OUT_FQ}"

echo "Done. Filtered FASTQ written to: ${OUT_FQ}"

