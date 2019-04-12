#create function to read, align, and calculate mean raw genetic distance (pi) from sequences.
library(stringr)
library(ape)
gen_dist_calc <- function(x, seq_type = "nexus") {
  if (seq_type == "nexus")
    rawseq <- read.nexus.data(x)
  else if (seq_type == "fasta")
    rawseq <- read.FASTA(x)
  binseq <- as.DNAbin(rawseq)
  alignseq <- clustalomega(binseq) 
  dist <- dist.dna(alignseq, model = "raw")
  cell <- str_split_fixed(x, pattern = "_", n = 4)[,4] #extract the cell id from the file name
  cell <- str_split_fixed(cell, pattern = "\\.", n = 2)[,1]
  cell <- str_split_fixed(cell, pattern = "_", n = 2)[,2]
  cell <- 
  list(cells =  cell, avg_pi = mean(dist), sd_pi = sd(dist))
}