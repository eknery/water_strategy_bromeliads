if(!require("ape")) install.packages("ape"); library("ape")

### input directory
dir_input = "0_data/"
### read acessions table
access = read.csv(paste0(dir_input, "taxon_access.csv") )

### loci names
all_loci = colnames(access)[colnames(access) != "taxon"]

### exrpoting directory
dir_out = "1_raw_sequences/"

for(locus_name in all_loci){
  ### lines with information
  access_clean = access[!is.na(access[,locus_name]),]
  ### accession numbers and names
  access_num = access_clean[,locus_name]
  access_names = access_clean[,"taxon"]
  ### downloading sequences
  one_locus = read.GenBank( access.nb = access_num )
  ### naming sequences
  names(one_locus) = access_names
  ### export
  write.dna(
    one_locus, 
    file = paste0(dir_out, locus_name, ".fasta"), 
    format = 'fasta' 
  )
  ### check
  print(paste0("Download complete: ", locus_name))
}

