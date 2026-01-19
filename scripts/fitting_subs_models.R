### load libraries
if(!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")

### input diretory
dir_input = "3_trimmed_sequences/"
file_names = list.files(dir_input)

### loading data
fasta_list = list()
for(i in 1:length(file_names) ){
  fasta_list[[i]] = read.phyDat(paste0(dir_input, file_names[i]),
                                format = "fasta",
                                type = "DNA"
  )
  names(fasta_list)[i] = str_remove(string = file_names[i], 
                                    pattern = ".fasta")
}

#### data frame for results
best_models = data.frame(
  locus = character(), 
  n = character(),
  seqLength = character(),
  bestModel = character(),
  stringsAsFactors = FALSE
)

### Loop through each alignment file
for (i in 1:length(fasta_list)) {
  ## one aligment
  one_fasta = fasta_list[[i]]
  ## fasta name
  locus = names(fasta_list)[i]
  ## number of sequences
  n = length(one_fasta)
  ## sequence length
  seqLength = length(one_fasta[[1]])
  ## fit models
  modelfits = modelTest(
    object = one_fasta, 
    model = "all", 
    G = TRUE, 
    k = 4, 
    I = TRUE
  )
  ## best fit model
  bestfit = modelfits[modelfits$AICcw == max(modelfits$AICcw),]
  ## Save result
  best_models = rbind(best_models, 
                   data.frame(locus = locus,
                              n = n,
                              seqLength = seqLength,
                              bestModel = bestfit$Model, 
                              stringsAsFactors = FALSE
                              )
                   )
}

### export dataframe
write.csv(best_models, 
          "4_model_fit/best_sub_models.csv", 
          row.names = FALSE
          )
