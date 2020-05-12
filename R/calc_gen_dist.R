#create function to read, align, and calculate mean raw genetic distance (pi) from fasta sequences.
library(stringr)
library(magrittr)
library(ape)
gen_dist_calc <- function(x) {
  alignseq <- read.FASTA(x)
  checkseq <- base.freq(alignseq, all = TRUE)
  dist <- dist.dna(alignseq, model = "raw", pairwise.deletion = TRUE) # this is Eq. 2 from Miraldo et al. 2016
  cell_1 <- str_split(x, pattern = "_", simplify = TRUE) #extract the cell id from the file name
  cell_2 <- cell_1[str_detect(cell_1, ".fas$")]
  cell <- str_split(cell_2, pattern = "\\.", simplify = TRUE)[,1]
  species_1 <- cell_1[str_detect(cell_1, "BOLD")]
  species <- str_split(species_1, pattern = ":", simplify = TRUE)[,2]
  list(species = species, 
       cell =  cell, 
       avg_pi = mean(dist), 
       sd_pi = sd(dist), 
       a = unname(checkseq["a"]),
       c = unname(checkseq["c"]),
       g = unname(checkseq["g"]),
       t = unname(checkseq["t"]),
       perc_missing = unname(checkseq["n"]),
       perc_gap = unname(checkseq["-"]),
       perc_ambiguous = unname(checkseq["?"])
       )
}
