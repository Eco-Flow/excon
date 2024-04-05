This is a Nextflow pipeline to run non-model organism Gene Ontology and CAFE 

When working with Human and Mouse, there are plentiful resources to take advantage of. Yet, for non-model organisms, it can take a long time to build the resources needed to analyse a new species. 

This pipeline aims to create a simple procedure to build a Gene ontology database for a new species. Then you can use it to produce GO enrichment hits with a few different programs to compare and contrast them.

It also optionally runs CAFE, so you can explore gene expansions and contractions


1. Running on mammalian data (default settings):

```
nextflow run main.nf -bg
```

Default settings:
params.ensembl_repo="genes"
params.ensembl_host='https://ensembl.org'
params.ensembl_dataset="example.txt"
params.focal = "Branchiostoma_lanceolatum.BraLan2.pep.all.fa"
params.predownloaded_fasta= "./Background_species_folder/*"
params.predownloaded_gofiles= "./Background_gofiles_folder/*"
params.outdir = "results"
params.download= false

To run on UCL myriad, you use the following command:
`nextflow run main.nf -bg -profile myriad -resume --focal Branchiostoma_lanceolatum.BraLan2.pep.all.fa`

`-resume`: ensure you continue from the last workable part of the pipeline.
`--focal Branchiostoma_lanceolatum.BraLan2.pep.all.fa`: explicitely says which file is the focal sample.

2. To test on an example (full) dataset using insect downloads:

```
nextflow run main.nf -profile myriad --ensembl_dataset example_insects.txt --download true -bg --focal My_polybia_prots --ensembl_repo metazoa_mart --ensembl_host https://metazoa.ensembl.org
```

#Notice you need to redirect to the insect ensembl section: with repo and host.
#You can work out which repo to use by following the later code below (section 2). Where you can search for samples on the ensembl database.
#Sometimes Ensembl biomart is down with error code 500, meaning you may need to download samples manually
`--download true`: tells the program to download the samples from ensembl






#Testing other features (in development)

1. To test orthofinder options use:
```
docker run -it --rm davidemms/orthofinder:2.5.4 orthofinder -h
docker run -it --rm chriswyatt/goatee_biomart R
```

To get a subset for testing this repo:
```
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pig.fasta | awk -F '\t' '{if(index($1,"Olfactory")!=0) printf("%s\n%s\n",$1,$2);}' | gzip > Pig_olfactory.fasta.gz
```

These lines are helpful to explore biomaRt :


2. To check what species are present and their species name codes:

```
library(biomaRt)
ensembl <- useEnsembl(biomart = "metazoa_mart", host="https://metazoa.ensembl.org")
datasets <- listDatasets(ensembl)
datasets
```

For mammals:
```
ensembl <- useEnsembl(biomart = "genes", host="https://ensembl.org")
```

Useful yml for future projects: https://github.com/nf-core/circrna/blob/e36d85792a9f9c2fc317ead0131560fddbae9462/environment.yml


Common errors:
```
Timeout was reached: [www.ensembl.org:80] Operation timed out after 300000 milliseconds with 11844514 bytes received
```
This error is due to biomaRt timing out. A way to counter this is to try the script again (with -resume), and to reduce the number of gene pulls in each loop (at line 9 of bin/R_biomart.R). From the default 200 to 100 or 50.
