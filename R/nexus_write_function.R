#function to write nexus files
nexus_write_fun <- function(x, out_folder) {
  for (i in seq_along(x)) {
    df <- x[[i]]
    for (j in seq_along(df)) {
      nuc.df <- as.data.frame(df[[j]]) #loop through each species per cell
      nuc.df <- nuc.df[order(sapply(nuc.df[,"nucleotides"], length), decreasing = TRUE),] #order df by nucleotide length (need to do this for padding short sequences with "N"s)
      nuc.vec <- strsplit(nuc.df[,"nucleotides"], "") #nucleotide column. Need to split nucleotides into individual characters.
      names(nuc.vec) <- paste(nuc.df[,"recordID"]) #names are the recordIDs
      nuc.pad <- pad_seqs_fun(nuc.vec) #pad nucleotides with "N"s
      write.nexus.data(nuc.pad, file = paste0(out_folder, "/", unique(nuc.df[,"bin_uri"]), "_", unique(nuc.df[,"cells"]), ".nex"), missing = "N", interleaved = FALSE) #write to a nexus file, which is named by the species and cell number
    }
  }
}