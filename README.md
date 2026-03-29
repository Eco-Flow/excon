# EXCON (v2.1.0)

A Nextflow pipeline to perform gene expansion and contraction analysis using CAFE.
It runs orthofinder to build input for CAFE (expansion and contraction analysis), then runs eggnogmapper to get GO terms, and run GO analysis.
Optionally, it also runs GO analysis on chromosomes, to search for enriched terms by chromosome.

It works with any set of species that have a genome (fasta) and annotation (gff) file. 
(minimum of 5 species ideally up to around 30). Maximum 100 species (normally). 

To run the GO annotation. You must provide a database yourself with `--eggnog_data_dir`, else everytime you run the pipeline, it will download the DB 
for you. So be careful, it is ~45GB. Please run it once, and save the DB somewhere handy to point to. 
This is then used to check what GO terms are associated with expanded or contracted gene sets (from CAFE).

## Overview

The general pipeline logic is as follows:

* Downloads genome and annotation files from NCBI `[NCBIGENOMEDOWNLOAD]`, or you provide your own.
* Unzips the files, if necessary `[GUNZIP]`
* Standardises and filters GFF annotations `[AGAT_CONVERTSPGXF2GXF]`.
* Extracts longest protein `[AGAT_SPKEEPLONGESTISOFORM]`.
* Gets the protein sequences `[GFFREAD]`.
* Renames the genes to gene name (as some will be isoform name) `RENAME_FASTA`.
* Optionally describes genome assembly and annotation:
  - `[BUSCO_BUSCO]`: Completeness of the genome compared to expected gene set.
  - `[QUAST]`: Assembly contiguity statistics (N50 etc).
  - `[AGAT_SPSTATISTICS]`: Gene, exon, and intron statistics.
* Finds orthologous genes across species `[ORTHOFINDER_CAFE]`.
* Rescales species tree branch lengths for CAFE `[RESCALE_TREE]`.
* Runs gene family evolution analysis `[CAFE]` and plots results `[CAFE_PLOT]`.
* Optionally downloads the eggnogmapper database `[EGGNOG_DOWNLOAD]`.
* Optionally assigns GO terms to genes using `[EGGNOGMAPPER]`.
* Optionally plots GO enrichment for expanded/contracted gene families `[CAFE_GO]`.
* Optionally plots GO enrichment of genes by chromosome `[CHROMO_GO]`.
* Optionally summarizes GO enrichment by chromosome`[SUMMARIZE_CHROMO_GO]`.

## Installation

Nextflow pipelines require a few prerequisites. There is further documentation on the nf-core webpage [here](https://nf-co.re/docs/usage/installation), about how to install Nextflow.

### Prerequisites

- [Docker](https://docs.docker.com/engine/install/) or [Singularity](https://docs.sylabs.io/guides/3.11/admin-guide/installation.html).
- [Java](https://www.java.com/en/download/help/download_options.html) and [openJDK](https://openjdk.org/install/) >= 8 (**Please Note:** When installing Java versions are `1.VERSION` so `Java 8` is `Java 1.8`).
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) >= `v25.10.0`.
- When running nextflow with this pipeline, ideally run `NXF_VER=25.10.0` beforehand, to ensure functionality on this version.

### Install

To install the pipeline please use the following commands but replace VERSION with a [release](https://github.com/Eco-Flow/excon/releases).

`wget https://github.com/Eco-Flow/excon/archive/refs/tags/VERSION.tar.gz -O - | tar -xvf -`

or

`curl -L https://github.com/Eco-Flow/excon/archive/refs/tags/VERSION.tar.gz --output - | tar -xvf -`

This will produce a directory in the current directory called `excon-VERSION` which contains the pipeline.

## Inputs

### Required

* `--input /path/to/csv/file` - A singular csv file as input in one of the two formats stated below.

This csv can take 2 forms:
* A 2 field csv where each row is a unique species name followed by a Refseq genome reference ID (**NOT** a Genbank reference ID) i.e. `data/input_small-s3.csv`. The pipeline will download the relevant genome fasta file and annotation gff3 (or gff augustus) file.
* A 3 field csv where each row is a unique species name, followed by an absolute path to a genome fasta file, followed by an absolute path to an annotation gff3 (or gff augustus) file. Input can be gzipped (.gz) or not.

**Please Note:** The genome has to be chromosome level not contig level.

2 fields (Name,Refseq_ID):
```
Drosophila_yakuba,GCF_016746365.2
Drosophila_simulans,GCF_016746395.2
Drosophila_santomea,GCF_016746245.2
```

3 fields (Name,genome.fna,annotation.gff):
```
Drosophila_yakuba,data/Drosophila_yakuba/genome.fna.gz,data/Drosophila_yakuba/genomic.gff.gz
Drosophila_simulans,data/Drosophila_simulans/genome.fna.gz,data/Drosophila_simulans/genomic.gff.gz
Drosophila_santomea,data/Drosophila_santomea/genome.fna.gz,data/Drosophila_santomea/genomic.gff.gz
```

> **Note:** Genomes should be chromosome-level, not contig-level. RefSeq IDs must be used (not GenBank IDs).

## Parameters

### Core options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Path to input CSV file | Required |
| `--outdir` | Output directory | `results` |
| `--groups` | NCBI taxonomy group for genome download (e.g. `insects`, `bacteria`) | `insects` |
| `--help` | Display help message | `false` |
| `--custom_config` | Path to a custom Nextflow config file | `null` |

### Quality statistics (optional)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--stats` | Run BUSCO, QUAST and AGAT statistics on genomes | `null` |
| `--busco_lineage` | BUSCO lineage database (e.g. `insecta_odb10`) | `null` |
| `--busco_mode` | BUSCO mode (`genome`, `proteins`, `transcriptome`) | `null` |
| `--busco_lineages_path` | Path to local BUSCO lineage databases | `null` |
| `--busco_config` | Path to BUSCO config file | `null` |

### CAFE gene family evolution

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--skip_cafe` | Skip CAFE analysis | `null` |
| `--cafe_max_differential` | Maximum gene count differential for CAFE filtering on retry | `50` |
| `--tree_scale_factor` | Scale factor for rescaling species tree branch lengths | `1000` |

### GO annotation with EggNOG-mapper (optional)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--run_eggnog` | Run EggNOG-mapper GO annotation | `false` |
| `--eggnog_data_dir` | Path to pre-downloaded EggNOG database directory | `null` |
| `--eggnog_target_taxa` | Restrict annotations to orthologs from this taxon and its descendants (NCBI taxon ID) | `null` |
| `--eggnog_tax_scope` | Taxonomic scope for orthologous group assignment (e.g. `50557` for Insecta) | `null` |
| `--eggnog_evalue` | Maximum e-value threshold for sequence matches | `null` |
| `--eggnog_score` | Minimum bitscore threshold for matches | `null` |
| `--eggnog_pident` | Minimum percent identity (%) | `null` |
| `--eggnog_query_cover` | Minimum query coverage (%) | `null` |
| `--eggnog_subject_cover` | Minimum subject coverage (%) | `null` |

**Note:** The EggNOG database is ~45GB. If `--eggnog_data_dir` is not provided, the database will be downloaded automatically on each run. We strongly recommend downloading it once and reusing it:
Then pass `--eggnog_data_dir /path/to/eggnog_data` to the pipeline.

### GO enrichment analysis (optional, requires `--run_eggnog`)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--chromo_go` | Run GO enrichment analysis by chromosome | `null` |
| `--go_cutoff` | P-value cutoff for GO enrichment | `0.05` |
| `--go_type` | GO test type (e.g. `none`) | `none` |
| `--go_max_plot` | Maximum number of GO terms to plot | `10` |

### Resource limits

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--max_memory` | Maximum memory per job | `128.GB` |
| `--max_cpus` | Maximum CPUs per job | `16` |
| `--max_time` | Maximum runtime per job | `240.h` |

## Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers |
| `singularity` | Run with Singularity containers |
| `conda` | Run with Conda environments |
| `test_bacteria` | Test run with small bacterial genomes |
| `test_small` | Test run with small insect genomes |


## Profiles

This pipeline is designed to run in various modes that can be supplied as a comma separated list i.e. `-profile profile1,profile2`.

### Container Profiles

Please select one of the following profiles when running the pipeline.

* `docker` - This profile uses the container software Docker when running the pipeline. This container software requires root permissions so is used when running on cloud infrastructure or your local machine (depending on permissions). **Please Note:** You must have Docker installed to use this profile.
* `singularity` - This profile uses the container software Singularity when running the pipeline. This container software does not require root permissions so is used when running on on-premise HPCs or you local machine (depending on permissions). **Please Note:** You must have Singularity installed to use this profile.
* `apptainer` - This profile uses the container software Apptainer when running the pipeline. This container software does not require root permissions so is used when running on on-premise HPCs or you local machine (depending on permissions). **Please Note:** You must have Apptainer installed to use this profile.

### Optional Profiles

* `local` - This profile is used if you are running the pipeline on your local machine.
* `aws_batch` - This profile is used if you are running the pipeline on AWS utilising the AWS Batch functionality. **Please Note:** You must use the `Docker` profile with with AWS Batch.
* `test_small` - This profile is used if you want to test running the pipeline on your infrastructure, running from predownloaded go files. **Please Note:** You do not provide any input parameters if this profile is selected but you still provide a container profile.
* `test_biomart` - This profile is used if you want to test running the pipeline on your infrastructure, running from the biomart input. **Please Note:** You do not provide any input parameters if this profile is selected but you still provide a container profile.

## Custom Configuration

If you want to run this pipeline on your institute's on-premise HPC or specific cloud infrastructure then please contact us and we will help you build and test a custom config file. This config file will be published to our [configs repository](https://github.com/Eco-Flow/configs). 

## Running the Pipeline

**Please note:** The `-resume` flag uses previously cached successful runs of the pipeline.

1. Example run the full test example data:

```
NXF_VER=25.04.8 # Is latest it is tested on
nextflow run main.nf -resume -profile docker,test_small
```

*Settings in test_small:*
input = "input_small-s3.csv"

For the fastest run use: `nextflow run main.nf -resume -profile docker,test_bacteria`

2. To run on your own data (minimal run), cafe only. 

```
# NXF_VER=25.04.8  
nextflow run main.nf -resume -profile docker --input data/input_small-s3.csv
```

3. To run with cafe and GO analysis

```
# NXF_VER=25.04.8
nextflow run main.nf -resume -profile docker --input data/input_small-s3.csv --chromo_go --go_type bonferoni --stats --run_eggnog --eggnog_data_dir /path/to/eggnogdb 
```

## Citation

This pipeline is published on Workflowhub using the nf-core template. If you use this pipeline in you work, the following citations are essential:

excon:
Christopher Wyatt, Gene EXpansion and CONtraction analysis pipeline (2026) https:/doi.org/10.48546/WORKFLOWHUB.WORKFLOW.2141.1.

nf-core:
Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.

If you used any of these tools within the pipeline, you must also cite:

CAFE:
Fábio K Mendes, Dan Vanderpool, Ben Fulton, Matthew W Hahn, CAFE 5 models variation in evolutionary rates among gene families, 
Bioinformatics, 2020; btaa1022, https://doi.org/10.1093/bioinformatics/btaa1022

Orthofinder:
Emms, D.M. and Kelly, S. (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics. 
Genome Biology 20:238

AGAT:
Dainat J. 2022. Another Gtf/Gff Analysis Toolkit (AGAT): Resolve interoperability issues and accomplish more with your annotations. Plant and Animal Genome XXIX Conference. https://github.com/NBISweden/AGAT.

eggNOG-mapper (if used):
Carlos P Cantalapiedra, Ana Hernández-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas, eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale, Molecular Biology and Evolution, Volume 38, Issue 12, December 2021, Pages 5825–5829, https://doi.org/10.1093/molbev/msab293

A full list of tools and their versions are found in the `software_versions.yml` in the `results/pipeline_info`. So ensure to look here for any additional tools you need to cite.

For the --stats module, you would need to cite BUSCO, Quast, AGAT

## Contact Us

If you need any support do not hesitate to contact us at any of:

`c.wyatt [at] ucl.ac.uk` 

`ecoflow.ucl [at] gmail.com`
