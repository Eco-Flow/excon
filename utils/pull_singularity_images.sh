#!/usr/bin/env bash
# Download every Singularity/Apptainer image excon uses into a cache directory,
# under the file names Nextflow looks for, so the pipeline can then run where
# there is no internet access (e.g. HPC compute nodes).
#
# Run it once on a machine with internet access (e.g. a login node), from any
# directory:
#
#   utils/pull_singularity_images.sh /shared/path/excon_images
#
# then point the pipeline at that directory, e.g. with
#   export NXF_SINGULARITY_LIBRARYDIR=/shared/path/excon_images
# (or NXF_APPTAINER_LIBRARYDIR with -profile apptainer). See docs/offline.md.
#
# Options:
#   --dry-run        list the images and the file names, download nothing
#   --engine NAME    'singularity' or 'apptainer' (default: whichever is installed)
#   --registry HOST  registry for images named without one (default: quay.io,
#                    as set by the pipeline's singularity/apptainer profiles)
#
# The image list comes from `nextflow inspect`, which lists the container of
# every process in the pipeline, used or not, so one cache serves any options.

set -euo pipefail

usage() { sed -n '2,21p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

dry_run=false
engine=""
registry="quay.io"
cache_dir=""

while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run)  dry_run=true ;;
        --engine)   engine="$2"; shift ;;
        --registry) registry="$2"; shift ;;
        -h|--help)  usage 0 ;;
        -*)         echo "ERROR: unknown option $1" >&2; usage 1 ;;
        *)          [ -z "$cache_dir" ] || { echo "ERROR: give one cache directory" >&2; usage 1; }
                    cache_dir="$1" ;;
    esac
    shift
done
[ -n "$cache_dir" ] || { echo "ERROR: give the directory to download the images into" >&2; usage 1; }

if [ -z "$engine" ] && ! $dry_run; then
    if command -v singularity >/dev/null 2>&1; then engine=singularity
    elif command -v apptainer >/dev/null 2>&1; then engine=apptainer
    else echo "ERROR: neither singularity nor apptainer is on the PATH" >&2; exit 1
    fi
fi
command -v nextflow >/dev/null 2>&1 || { echo "ERROR: nextflow is not on the PATH" >&2; exit 1; }
command -v python3  >/dev/null 2>&1 || { echo "ERROR: python3 is not on the PATH" >&2; exit 1; }

pipeline_dir="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

echo "Listing the pipeline's containers with nextflow inspect ..." >&2
# inspect prints the JSON, then the pipeline's own completion line.
( cd "$work" && nextflow inspect "$pipeline_dir/main.nf" -profile "${engine:-singularity}" --outdir "$work/out" ) \
    > "$work/inspect.txt"

# Resolve each container as Nextflow does (ContainerHandler.normalizeSingularityImageName
# and SingularityCache.simpleName): a name with no scheme is a Docker image, prefixed
# with the registry unless its first path part is a host; the cached file name is the
# URL without its scheme, with ':' and '/' replaced by '-', plus '.img'.
python3 - "$work/inspect.txt" "$registry" > "$work/images.tsv" <<'PY'
import json, re, sys
text = open(sys.argv[1]).read()
data, _ = json.JSONDecoder().raw_decode(text[text.index('{'):])
registry = sys.argv[2].rstrip('/')
seen = set()
for proc in data["processes"]:
    img = proc.get("container")
    if not img or img in seen:
        continue
    seen.add(img)
    if not re.match(r'^[^/:. ]+://', img):
        first = img.split('/', 1)[0] if '/' in img else ''
        absolute = '/' in img and ('.' in first or ':' in first)
        img = "docker://" + (img if absolute or not registry else registry + "/" + img)
    name = img.split('://', 1)[1]
    ext = '.img'
    if '.sif:' in name:
        ext, name = '.sif', name.replace('.sif:', '-')
    elif name.endswith('.sif'):
        ext, name = '.sif', name[:-4]
    print(img + "\t" + name.replace(':', '-').replace('/', '-') + ext)
PY

n=$(wc -l < "$work/images.tsv" | tr -d ' ')
echo "$n images" >&2
mkdir -p "$cache_dir"

failed=0
while IFS=$'\t' read -r url file; do
    if [ -f "$cache_dir/$file" ]; then
        echo "have     $file"
    elif $dry_run; then
        echo "would pull $url -> $file"
    else
        echo "pulling  $url"
        # Pull to a temporary name first, so an interrupted download never
        # leaves a partial image that Nextflow would take for a complete one.
        if "$engine" pull --name "$cache_dir/$file.part" "$url"; then
            mv "$cache_dir/$file.part" "$cache_dir/$file"
        else
            rm -f "$cache_dir/$file.part"
            echo "FAILED   $url" >&2
            failed=$((failed + 1))
        fi
    fi
done < "$work/images.tsv"

if [ "$failed" -gt 0 ]; then
    echo "ERROR: $failed of $n images failed to download; re-run to retry them." >&2
    exit 1
fi
$dry_run || echo "All $n images are in $cache_dir" >&2
