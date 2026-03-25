# Make a dir and mv to it. This part user needs to change.
mkdir -p ~/Downloads/eggnog_data
cd ~/Downloads/eggnog_data

# Correct URLs using eggnog5.embl.de
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz

# Decompress
gunzip eggnog.db.gz
gunzip eggnog_proteins.dmnd.gz
tar -xzf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz