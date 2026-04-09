#!/usr/bin/env bash
# =============================================================================
# validate_inputs.sh
# Checks every GCF/GCA accession in all benchmark input CSVs against NCBI.
# Prints OK / FAIL for each, and a summary at the end.
#
# Requires: ncbi-genome-download (pip install ncbi-genome-download)
#           OR the NCBI datasets CLI (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
#
# Usage:
#   ./benchmark/validate_inputs.sh                  # check all CSVs
#   ./benchmark/validate_inputs.sh inputs/bacteria_close_contiguous.csv
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUTS_DIR="${SCRIPT_DIR}/inputs"

# Which CSVs to check
if [[ $# -gt 0 ]]; then
    CSV_FILES=("$@")
else
    CSV_FILES=("${INPUTS_DIR}"/*.csv)
fi

# Detect which tool is available
if command -v datasets &>/dev/null; then
    TOOL="datasets"
elif command -v ncbi-genome-download &>/dev/null; then
    TOOL="ngd"
else
    echo "ERROR: Neither 'datasets' (NCBI CLI) nor 'ncbi-genome-download' found."
    echo "Install one of:"
    echo "  conda install -c conda-forge ncbi-datasets-cli"
    echo "  pip install ncbi-genome-download"
    exit 1
fi

echo "Using: $TOOL"
echo ""

total=0
failed=0
failed_list=()

check_accession_datasets() {
    local acc="$1"
    # datasets summary returns JSON; non-zero exit or empty total means not found
    local result
    result=$(datasets summary genome accession "$acc" 2>/dev/null)
    local count
    count=$(echo "$result" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('total_count',0))" 2>/dev/null || echo 0)
    [[ "$count" -gt 0 ]]
}

check_accession_ngd() {
    local acc="$1"
    # --dry-run exits 0 if accession found, non-zero if not
    ncbi-genome-download --dry-run -A "$acc" all &>/dev/null
}

for csv in "${CSV_FILES[@]}"; do
    echo "=== $(basename "$csv") ==="
    while IFS=',' read -r species accession; do
        [[ -z "$species" || -z "$accession" ]] && continue
        accession="${accession//[$'\r\n']}"   # strip Windows line endings
        ((total++)) || true

        if [[ "$TOOL" == "datasets" ]]; then
            ok=$(check_accession_datasets "$accession" && echo yes || echo no)
        else
            ok=$(check_accession_ngd "$accession" && echo yes || echo no)
        fi

        if [[ "$ok" == "yes" ]]; then
            printf "  OK   %-45s %s\n" "$species" "$accession"
        else
            printf "  FAIL %-45s %s\n" "$species" "$accession"
            ((failed++)) || true
            failed_list+=("$(basename "$csv"): ${species},${accession}")
        fi
    done < "$csv"
    echo ""
done

echo "════════════════════════════════════════"
echo " Total: ${total}  OK: $((total - failed))  FAILED: ${failed}"
echo "════════════════════════════════════════"

if [[ ${#failed_list[@]} -gt 0 ]]; then
    echo ""
    echo "Failed accessions to fix:"
    for entry in "${failed_list[@]}"; do
        echo "  $entry"
    done
    exit 1
fi
