process EGGNOG_TO_OG_GO {
    tag "og_go"
    label 'process_single'
    label 'process_long'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
        path(go_files)
        path(n0_tsv)

    output:
        path("OG_GO_format.tsv"),                                                                        emit: og_go
        tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    script:
    """
    python3 <<'EOF'
import glob, os, sys

print("Reading N0.tsv...", flush=True)
with open("N0.tsv") as fh:
    header = fh.readline().strip().split("\\t")

    fixed_cols = {"HOG", "OG", "Gene Tree Parent Clade", "Orthogroup"}
    species_start = next(i for i, h in enumerate(header) if h.split(".")[0] not in fixed_cols)
    species_cols = [h.split(".")[0] for h in header[species_start:]]
    print(f"Found {len(species_cols)} species columns", flush=True)

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
    print(f"Found {len(og_genes)} orthogroups", flush=True)

needed_genes = set()
for genes in og_genes.values():
    needed_genes.update(genes)
print(f"Found {len(needed_genes)} unique (species, gene) pairs needed", flush=True)

gene_go = {}
go_files = glob.glob("*.go.txt")
print(f"Reading {len(go_files)} GO files...", flush=True)
for i, go_file in enumerate(go_files):
    sp = go_file.replace(".go.txt", "")
    with open(go_file) as fh:
        for line in fh:
            parts = line.strip().split("\\t")
            if len(parts) == 2:
                gene, go = parts
                if (sp, gene) in needed_genes:
                    gene_go.setdefault((sp, gene), set()).add(go)
    print(f"  [{i+1}/{len(go_files)}] {sp}", flush=True)

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
