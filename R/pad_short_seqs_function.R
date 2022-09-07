#####Code for padding short sequences with gaps
#function to pad the ends of shorter sequences with "N"s so all sequences are the same length for alignments
pad_seqs_fun <- function(nuc) {
  nuc <- nuc[order(sapply(nuc, length), decreasing = TRUE)] #order sequences in descending order by length. This removes the necessity to do a bunch of unnecessary pairwise comparisons among sequences. Hashing this out because I'm ordering the data frame by sequence length before I apply this function. 
  nuc.out <- vector("list", length(nuc)) #establish an empty vector to fill
  nuc.out[[1]] <- nuc[[1]] #replace the first sequence in the out vector
  names(nuc.out) <- names(nuc) #keep the names of the input vector
  for (i in seq(2, length(nuc)))
    if(length(nuc[[i]]) <= length(nuc[[1]])) #if the sequence is shorter than the longest sequence, append the missing data designator "N" to the end of the sequence
      nuc.out[[i]] <- append(nuc[[i]], rep(c("N"), length(nuc[[1]]) - length(nuc[[i]])))
  return(nuc.out)
}