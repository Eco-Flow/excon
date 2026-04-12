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
    # datasets summary returns JSON — check accession exists AND has annotation
    local result
    result=$(datasets summary genome accession "$acc" 2>/dev/null)
    local count has_annotation
    count=$(echo "$result" | python3 -c "
import sys, json
d = json.load(sys.stdin)
print(d.get('total_count', 0))
" 2>/dev/null || echo 0)
    [[ "$count" -gt 0 ]] || return 1

    # Check annotation_info is present (means GFF exists on NCBI)
    has_annotation=$(echo "$result" | python3 -c "
import sys, json
d = json.load(sys.stdin)
reports = d.get('reports', [])
has = any(r.get('annotation_info') for r in reports)
print('yes' if has else 'no')
" 2>/dev/null || echo "unknown")

    if [[ "$has_annotation" == "no" ]]; then
        # Annotate the failure reason for the caller
        FAIL_REASON="no annotation (GFF) on NCBI"
        return 1
    fi
    return 0
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

        FAIL_REASON=""
        if [[ "$TOOL" == "datasets" ]]; then
            ok=$(check_accession_datasets "$accession" && echo yes || echo no)
        else
            ok=$(check_accession_ngd "$accession" && echo yes || echo no)
        fi

        if [[ "$ok" == "yes" ]]; then
            printf "  OK   %-45s %s\n" "$species" "$accession"
        else
            reason="${FAIL_REASON:-not found on NCBI}"
            printf "  FAIL %-45s %s  [%s]\n" "$species" "$accession" "$reason"
            ((failed++)) || true
            failed_list+=("$(basename "$csv"): ${species},${accession}  [${reason}]")
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
