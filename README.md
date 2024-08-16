A Nextflow pipeline to run gene expansion and contraction analysis with CAFE. 

It works with any set of species that have a genome (fasta) and annotation (gff) file. 
(minimum of 5 species ideally up to around 15).

As well as running CAFE, you can also run GO enrichment analysis (with user-provided GO files, or with GO files semi-automatically downloaded from Ensembl biomart). This will check what GO terms are associated with expanded or contracted gene sets.


The general pipeline logic is as follows:

* Downloads the genome and gene annotation files from NCBI `[DOWNLOAD]`.
  - Or you provide your own genomes/annotations
* Extract longest protein fasta sequences `[GFFREAD]`.
* Finds orthologous genes `[ORTHOFINDER_CAFE]`.
* Runs cafe analysis on the output of orthofinder `[CAFE]`.
* Runs gene to GO assignment (optional) `[ORTHOFINDER_GO], [GO_ASSIGN]`.  REQUIRES GO FILES TO RUN, check [optional inputs](#optional).
* Plot GO enrichment for excon genes (optional) `[CAFE_GO]`.

## Installation

Nextflow pipelines require a few prerequisites. There is further documentation on the nf-core webpage [here](https://nf-co.re/docs/usage/installation), about how to install Nextflow.

### Prerequisites

- [Docker](https://docs.docker.com/engine/install/) or [Singularity](https://docs.sylabs.io/guides/3.11/admin-guide/installation.html).
- [Java](https://www.java.com/en/download/help/download_options.html) and [openJDK](https://openjdk.org/install/) >= 8 (**Please Note:** When installing Java versions are `1.VERSION` so `Java 8` is `Java 1.8`).
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) >= `v23.07.0`.

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

### Optional 

* `--outdir /path/to/output/directory` - A path to the output directory where the results will be written to (**Default:** `Results`).
* `--predownloaded_fasta /path/to/predownloaded/fasta` - A path to a folder containing a singular or multiple file for each species you wish to use for GO assignment (**Default:** `null`).
* `--predownloaded_gofiles /path/to/predownloaded_gofiles` - A path to a folder containing a singular or multiple file for each species you wish to use for GO assignment. Matched to above flag. One GO file and one protein fasta for each species with matched names (**Default:** `null`). 
* `--ensembl_dataset  /path/to/ensembl/dataset` - A path to a file containing a list of ensembl biomart dataset names separated by newline i.e. `data/biomart.txt` (**Default:** `null`).
* `--ensembl_biomart /path/to/ensembl_biomart` - A name of an ensembl biomart that contains the species data you want (e.g. metazoa_mart for insects) (**Default:** `null`).
* `--skip_cafe` - A flag to skip the cafe section. Used if you just wish to run go assignment for a species without runnig Cafe (**Default:** `null`).
* `--go_expansion` - A flag optionally choose to run a basic expansion/contraction analysis (**Default:** `null`).
* `--chromo_go` - A flag to optionally choose to run GO analysis on each chromosome (**Default:** `null`).
* `--clean` - A true or false value assigned to this parameter will determine whether the work directory is automatically deleted or not if the pipeline is successful. Deleting the work directory saves space however it will not be possible to use this work directory in future for caching (**Default:** `false`).
* `--help` - A true value assigned to this parameter will cause the help message to be displayed instead of pipeline running (**Default:** `false`).
* `--custom_config` - A path or URL to a custom configuration file.
* `--busco` - A flag to optionally choose to run BUSCO on each genome.
* `--atag` - A flag to optionally choose to run ATAG statistics on each genome.

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
nextflow run main.nf -resume -profile docker,test_small
```

*Settings in test_small:*
input = "input_small-s3.csv"
predownloaded_fasta = "s3://excon/data/Insect_data/fasta/*"
predownloaded_gofiles = "s3://excon/data/Insect_data/gofiles/*"

2. To run on your own data (minimal run), cafe only. 

```
nextflow run main.nf -resume -profile docker --input data/input_small-s3.csv
```

3. To run on your own data with GO enrichment analysis (using predownloaded fasta/go files for GO assignment)

```
nextflow run main.nf -resume -profile docker --input data/input_small-s3.csv --cafe \|
--predownloaded_fasta 's3://excon/data/Insect_data/fasta/*' --predownloaded_gofiles 's3://excon/data/Insect_data/gofiles/*' 
```

4. To run on your own data with GO enrichment analysis + retrieval of GO assignment species

If you do not have GO files to run GO enrichment, you can run the following code to semi-auto download them from NCBI biomart.

You first need to go to Ensembl Biomart to find the species IDs you want to use to assign GO terms to your species. Ideally you should choose one or more species that are closely related and have good GO annotations.

i) To check what species are present and their species name codes you need to download the biomaRt library in R (for metazoa):

```
library(biomaRt)
ensembl <- useEnsembl(biomart = "metazoa_mart", host="https://metazoa.ensembl.org")
datasets <- listDatasets(ensembl)
datasets
```

You will see something like:

```
                       dataset
1     aagca019059575v1_eg_gene
2       aagca914969975_eg_gene
3     aagca933228735v1_eg_gene
4           aalbimanus_eg_gene
5          aalbopictus_eg_gene
6            aalvpagwg_eg_gene
```

The dataset IDs are what you need to enter into the Nextflow script.

For mammals:
```
ensembl <- useEnsembl(biomart = "genes", host="https://ensembl.org")
```

Then you can run the excon script as follows:

```
nextflow run main.nf -resume -profile <apptainer/docker/singularity> --input data/input_small-s3.csv --cafe --cafe_go --ensembl_biomart "metazoa_mart" --ensembl_dataset "example.txt"
```

where `example.txt` is a file of dataset IDs from ensembl biomart (as shown above), separated by newline characters.

e.g.:

```
aagca019059575v1_eg_gene
aagca914969975_eg_gene
aagca933228735v1_eg_gene
```

## Citation

This pipeline is not yet published. Please contact us if you wish to use our pipeline, we are happy to help you run the pipeline.

## Contact Us

If you need any support do not hesitate to contact us at any of:

`c.wyatt [at] ucl.ac.uk` 

`ecoflow.ucl [at] gmail.com`
