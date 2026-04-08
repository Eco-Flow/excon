process OG_ANNOTATION_SUMMARY {
    tag "og_annotation_summary"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    path annotations   // collected *.emapper.annotations files (one per species)
    path orthogroups   // N0.tsv or Orthogroups.tsv from OrthoFinder
    val  rep_species   // "" = auto-select; otherwise the species id to use

    output:
    path "OG_annotation_summary.tsv", emit: summary
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    script:
    """
    python3 <<'EOF'
import glob, os

# ---------------------------------------------------------------------------
# 1. Parse all annotation files into a dict: species -> {gene -> row_dict}
# ---------------------------------------------------------------------------
# Annotation columns (0-based):
#   0  query        1  seed_ortholog   2  evalue   3  score
#   4  eggNOG_OGs   5  max_annot_lvl  6  COG_cat  7  Description
#   8  Preferred_name  9  GOs  10  EC  11  KEGG_ko  ...  last  PFAMs

def parse_annot_file(path):
    sp = os.path.basename(path).replace(".emapper.annotations", "")
    genes = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 9:
                continue
            gene      = cols[0]
            score     = float(cols[3]) if cols[3] not in ("-", "") else 0.0
            cog       = cols[6]  if cols[6]  not in ("-", "") else ""
            desc      = cols[7]  if cols[7]  not in ("-", "") else ""
            pref      = cols[8]  if cols[8]  not in ("-", "") else ""
            kegg_ko   = cols[11] if len(cols) > 11 and cols[11] not in ("-", "") else ""
            pfams     = cols[-1] if cols[-1] not in ("-", "") else ""
            genes[gene] = dict(score=score, cog=cog, desc=desc, pref=pref,
                               kegg_ko=kegg_ko, pfams=pfams)
    return sp, genes

all_sp_annots = {}   # species -> {gene -> info}
for f in sorted(glob.glob("*.emapper.annotations")):
    sp, genes = parse_annot_file(f)
    all_sp_annots[sp] = genes

print(f"Loaded annotations for {len(all_sp_annots)} species", flush=True)

# ---------------------------------------------------------------------------
# 2. Pick the representative species
# ---------------------------------------------------------------------------
forced = "${rep_species}".strip()

if forced:
    rep = forced
    if rep not in all_sp_annots:
        raise ValueError(
            f"--eggnog_rep_species '{rep}' not found. "
            f"Available: {sorted(all_sp_annots.keys())}"
        )
    print(f"Using forced representative species: {rep}", flush=True)
else:
    # Auto: species with the most genes that have a non-empty description
    rep = max(all_sp_annots,
              key=lambda s: sum(1 for g in all_sp_annots[s].values() if g["desc"]))
    print(f"Auto-selected representative species: {rep} "
          f"({sum(1 for g in all_sp_annots[rep].values() if g['desc'])} annotated genes)",
          flush=True)

rep_annots = all_sp_annots[rep]

# ---------------------------------------------------------------------------
# 3. Read orthogroups table
# ---------------------------------------------------------------------------
print("Reading orthogroups table...", flush=True)
with open("${orthogroups}") as fh:
    header = fh.readline().strip().split("\\t")

fixed_cols = {"HOG", "OG", "Gene Tree Parent Clade", "Orthogroup"}
species_start = next(
    (i for i, h in enumerate(header) if h.split(".")[0] not in fixed_cols),
    1
)
species_cols = [h.split(".")[0] for h in header[species_start:]]

# Find column index for the representative species
rep_col_idx = None
for i, sp in enumerate(species_cols):
    if sp == rep:
        rep_col_idx = i
        break

if rep_col_idx is None:
    raise ValueError(
        f"Representative species '{rep}' not found as a column in the orthogroups table. "
        f"Columns: {species_cols[:5]}..."
    )

print(f"  Representative species column index: {rep_col_idx}", flush=True)

# ---------------------------------------------------------------------------
# 4. For each OG, find the best-scoring gene from the rep species
# ---------------------------------------------------------------------------
print("Writing OG_annotation_summary.tsv...", flush=True)

with open("${orthogroups}") as fh, open("OG_annotation_summary.tsv", "w") as out:
    fh.readline()  # skip header
    out.write("HOG\\trep_species\\trep_gene\\tdescription\\tpreferred_name\\tCOG_category\\tKEGG_ko\\tPFAMs\\n")

    for line in fh:
        cols = line.rstrip("\\n").split("\\t")
        og   = cols[0]

        # Genes from the representative species in this OG
        genes_str = cols[species_start + rep_col_idx] if (species_start + rep_col_idx) < len(cols) else ""
        genes = [g.strip() for g in genes_str.split(",") if g.strip()]

        # Pick best-scoring annotated gene; fall back to any gene present
        best_gene = None
        best_info = None
        best_score = -1.0

        for gene in genes:
            info = rep_annots.get(gene)
            if info and info["score"] > best_score:
                best_score = info["score"]
                best_gene  = gene
                best_info  = info

        if best_gene is None and genes:
            # Gene present in OG but not annotated
            best_gene  = genes[0]
            best_info  = dict(cog="", desc="", pref="", kegg_ko="", pfams="")

        if best_gene is None:
            # Rep species has no genes in this OG
            out.write(f"{og}\\t{rep}\\t-\\t-\\t-\\t-\\t-\\t-\\n")
        else:
            out.write(
                f"{og}\\t{rep}\\t{best_gene}\\t"
                f"{best_info['desc']}\\t{best_info['pref']}\\t"
                f"{best_info['cog']}\\t{best_info['kegg_ko']}\\t"
                f"{best_info['pfams']}\\n"
            )

print("Done.", flush=True)
EOF
    """

    stub:
    """
    echo -e "HOG\\trep_species\\trep_gene\\tdescription\\tpreferred_name\\tCOG_category\\tKEGG_ko\\tPFAMs" > OG_annotation_summary.tsv
    """
}
