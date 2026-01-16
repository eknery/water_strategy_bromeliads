
data = read.table(
  "0_data/data.csv",
  sep=",",
  h=T)

access = read.table(
  "0_data/genbank_accessions.csv",
  sep=",",
  h=T)

############################### DATA PROCESSING ################################


## genera and their identifiers
### split string
sp_str = strsplit(access$species, split = " ")
### get identifiers column
taxon = c()
for(i in 1:length(sp_str)){
  x = paste0(sp_str[[i]][1], "_", sp_str[[i]][2])
  taxon = c(taxon, x)
}
### adding identifier
access$taxon = taxon