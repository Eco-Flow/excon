library(biomaRt)

#List all deafult ensembl mart dbs
listEnsembl()
# make a db with the "genes" biomart
ensembl <- useEnsembl(biomart = "genes")
# Check out what is in this db
datasets <- listDatasets(ensembl)
head(datasets)

#Search for a term within the open DB
searchDatasets(mart = ensembl, pattern = "hsapiens")
# Select a dataset to work with.
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
#Choose a dataset in one go.
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#Use a different version of ensembl (could be useful for insects)
listEnsemblArchives()
# Then try a different version:
ensembl <- useEnsembl(biomart = "genes", version=80)

#Check out other species genome sites:
listEnsemblGenomes()

#Then find a specific species , e.g. APis:
ensembl_metazoa <- useEnsemblGenomes(biomart = "metazoa_mart")
searchDatasets(ensembl_metazoa, pattern = "Apis")

#Then select apis:
ensembl_apis <- useEnsemblGenomes(biomart = "metazoa_mart", dataset="amellifera_eg_gene")

#Now we can explore the filters present for the account:
filters = listFilters(ensembl_apis)
attributes = listAttributes(ensembl_apis)


getBM(attributes = c('with_ensembl_metazoa', 'with_go'), mart = "ensembl_apis")

#Doesn't work for insects , gruffffffff!!!

#But in human:
getBM(attributes = c('entrezgene_id','hgnc_symbol'), 
      filters = 'go', 
      values = 'GO:0004707', 
      mart = ensembl)

cannon<- getBM(attributes = c('ensembl_gene_id','go_id'), mart = ensembl)

# the above timed out so did (with a mirror, which now works):
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror="uswest")
cannon<- getBM(attributes = c('ensembl_gene_id','go_id'), mart = ensembl)

#Try again apis:
ensembl_apis <- useEnsembl(biomart = "metazoa_mart", dataset="amellifera_eg_gene", mirror="uswest" , host="https://metazoa.ensembl.org")
bazooka <- getBM(attributes = c('ensembl_gene_id','go_id'), mart = ensembl_apis)
write.table(bazooka, "go_hash.txt", row.names=T, quote=F, sep="\t" )

bazooka2 <- getBM(attributes = c('ensembl_gene_id','peptide'), mart = ensembl_apis)
#Issue now is that the peptide one times out before the job is done, a way around this is to provide a list of ids, and then request the sequecnes separately, in a loop before joining everything together :
affyids <- c("LOC724736")
bazooka2 <- getBM(attributes = c('ensembl_gene_id','peptide'), mart = ensembl_apis, values=affyids, filters='ensembl_gene_id',  uniqueRows=T)


#new script to split up the seqs:
bazooka <- getBM(attributes = c('ensembl_gene_id','go_id'), mart = ensembl_apis)
my_ids<-unique(bazooka$ensembl_gene_id)
d<-my_ids
new_ids<- split(d, ceiling(seq_along(d)/1000))

for(i in 1:length(new_ids)) { res_[i] <- getBM(attributes = c('ensembl_gene_id','peptide'), mart = ensembl_apis, values= new_ids[i], filters='ensembl_gene_id',  uniqueRows=T) }

#try
cabbage=new_ids[1]
for(i in 1:length(new_ids)) { res_[i] <- getBM(attributes = c('ensembl_gene_id','peptide'), mart = ensembl_apis, values= inputs, filters='ensembl_gene_id',  uniqueRows=T) }

new_ids2[i]$`10`

#works

for(i in 1:4) {
+ bazooka2 <- getBM(attributes = c('ensembl_gene_id','peptide'), mart = ensembl_apis, values=new_ids2[i], filters='ensembl_gene_id',  uniqueRows=T)
+ outname <- paste("Myoutput_", i, ".fasta")
+ write.table(bazooka2, )
append=        col.names=     dec=           eol=           file=          fileEncoding=  na=            qmethod=       quote=         row.names=     sep=           x=
+ write.table(bazooka2, file=outname)
+ }