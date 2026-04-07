process EGGNOG_TO_OG_GO {
    tag "og_go"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    path go_files      // per-species *.go.txt files
    path orthogroups   // N0.tsv or Orthogroups.tsv from OrthoFinder

    output:
        path("OG_GO_format.tsv"),                                                                        emit: og_go
        tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    script:
    """
    python3 <<'EOF'
import glob, os, sys
from concurrent.futures import ProcessPoolExecutor, as_completed

CPUS = ${task.cpus}

# ---------------------------------------------------------------------------
# 1. Read orthogroups table (N0.tsv or Orthogroups.tsv)
# ---------------------------------------------------------------------------
print("Reading orthogroups table...", flush=True)
with open("${orthogroups}") as fh:
    header = fh.readline().strip().split("\\t")

    fixed_cols = {"HOG", "OG", "Gene Tree Parent Clade", "Orthogroup"}
    species_start = next(
        (i for i, h in enumerate(header) if h.split(".")[0] not in fixed_cols),
        1  # Orthogroups.tsv has no fixed prefix cols; data starts at col 1
    )
    species_cols = [h.split(".")[0] for h in header[species_start:]]
    print(f"  {len(species_cols)} species columns", flush=True)

    og_genes = {}
    for line in fh:
        parts = line.strip().split("\\t")
        og = parts[0]
        gene_list = []
        for sp, genes_str in zip(species_cols, parts[species_start:]):
            for gene in genes_str.split(", "):
                gene = gene.strip()
                if gene:
                    gene_list.append((sp, gene))
        og_genes[og] = gene_list

print(f"  {len(og_genes)} orthogroups", flush=True)

needed_genes = set()
for genes in og_genes.values():
    needed_genes.update(genes)
print(f"  {len(needed_genes)} unique (species, gene) pairs needed", flush=True)

# ---------------------------------------------------------------------------
# 2. Read GO files in parallel
# ---------------------------------------------------------------------------
def read_go_file(go_file):
    sp = os.path.basename(go_file).replace(".go.txt", "")
    result = {}
    with open(go_file) as fh:
        for line in fh:
            parts = line.strip().split("\\t")
            if len(parts) == 2:
                gene, go_term = parts
                key = (sp, gene)
                if key in needed_genes:
                    result.setdefault(key, set()).add(go_term)
    return result

go_files = glob.glob("*.go.txt")
print(f"Reading {len(go_files)} GO files with {CPUS} workers...", flush=True)

gene_go = {}
with ProcessPoolExecutor(max_workers=CPUS) as pool:
    futures = {pool.submit(read_go_file, f): f for f in go_files}
    for i, fut in enumerate(as_completed(futures), 1):
        for key, gos in fut.result().items():
            gene_go.setdefault(key, set()).update(gos)
        print(f"  [{i}/{len(go_files)}] {os.path.basename(futures[fut])}", flush=True)

# ---------------------------------------------------------------------------
# 3. Write OG_GO_format.tsv
# ---------------------------------------------------------------------------
print("Writing OG_GO_format.tsv...", flush=True)
with open("OG_GO_format.tsv", "w") as out:
    for og, genes in og_genes.items():
        gos = set()
        for sp_gene in genes:
            gos.update(gene_go.get(sp_gene, set()))
        for go in sorted(gos):
            out.write(f"{og}\\t{go}\\n")

print("Done.", flush=True)
EOF
    """
}
