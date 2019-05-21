#create function to read, align, and calculate mean raw genetic distance (pi) from fasta sequences.
library(stringr)
library(magrittr)
library(ape)
gen_dist_calc <- function(x) {
  alignseq <- read.FASTA(x)
  dist <- dist.dna(alignseq, model = "raw")
  cell <- str_split_fixed(x, pattern = "_", n = 5)[,5] #extract the cell id from the file name
  cell <- str_split_fixed(cell, pattern = "\\.", n = 2)[,1]
  species <- str_split_fixed(cell, pattern = "_", n = 2)[,1]
  species <- str_split_fixed(species, pattern = ":", n = 2)[,2]
  cell <- str_split_fixed(cell, pattern = "_", n = 2)[,2]
  list(species = species, cell =  cell, avg_pi = mean(dist), sd_pi = sd(dist))
}