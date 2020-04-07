#' ---
#' title: "BigAss BOLD"
#' always_allow_html: yes
#' output: # use rmarkdown::render(<your-rmd-file.rmd>, output_format ="all" in the console to render both outputs
#'   html_notebook:
#'     theme: flatly
#'     highlight: tango
#'   github_document:
#'     toc: true
#' ---
#' Load packages, read in data, and filter for NAs
## ----include = FALSE-----------------------------------------------------
library(data.table)
#library(furrr)
library(leaflet)
library(sf)
library(raster)
library(ape)
library(gdata)
library(entropy)
library(here)
library(tidyverse)


#function to pad the ends of shorter sequences with "N"s so all sequences are the same length for alignments
source("R/pad_short_seqs_function.R") #pad short sequences with Ns
source("R/fasta_write_function.R") #write fasta files 
source("R/calc_gen_dist.R") #align seqs and calculate mean pi per species per cell
source("R/hill_num_calc.R") #calculate hill number of pi. Note* using clustal omega in the ape R package requires having a copy of clustal omega downloaded and accessible either through their PATH or have it indicated in the clustalomega() argument.
source("R/pi_summary_function.R")
source("R/sumstat_plot_function.R")

#### Make output folders ----
#will throw warnings if you're completing steps 1 and 2 on the same day. You can ignore these.

#if you want to work out of a folder from an earlier date, replace this string with the date
todays_date <- Sys.Date()

#folder for the entire project's output to go into
todays_results <- paste0("results_", todays_date)
dir.create(todays_results)

#output folder for spreadsheets, tables, etc.
output <- paste0(todays_results, "/output")
dir.create(output)

#folder for plots
plots <- paste0(todays_results, "/plots")
#make plots folder
dir.create(plots)

#folder for rasters
rasters <- paste0(todays_results, "/rasters")
#make raster folder
dir.create(rasters)
##### Align and calculate summary statistics on the sequences.

#read in the test nuc data set. If you're redoing this after filtering the data on a different day, replace the path with the correct test_nuc_*.csv path
test_nuc_path <- paste0(todays_results,"/output/test_nuc.csv")
test_nuc <- fread(test_nuc_path)


#split data frame by cell number and species
species_seq_split_one <- test_nuc %>%
  as.data.frame() %>% #some functions don't like spatial data frames
  drop.levels() %>%
  split(.$cells) %>% #split into a list of data frames, grouped by cell. Splitting by both cells and species at once doesn't work.
  lapply(drop.levels) %>% #drop any levels in the data frame. Have to perform first because extra factor levels can mess up the split function
  lapply(function(x){
    split(x, x$bin_uri)}) %>% #split each cell by species
  drop.levels()

#create directory to put fasta files
fasta_folder <- paste0(todays_results, "/output/bold_seqs")
dir.create(fasta_folder)

#write fasta files to output folder
fasta_write_fun(species_seq_split_one, out_folder = fasta_folder)

print("Wrote fasta files")

#' 
#' Align sequences. Have to execute this outside of the R environment
#' This takes way too long to do in R
#' I couldn't figure out an easy way to parallelize the process, so I came up with an ugly workaround
#' In the other_scripts folder I have two scripts, split_folders.txt and alignment_bash_script.txt
#' The split_folders.txt file contains a script to split the unaligned sequences into multiple folders. 
#' You can adjust the number of folders by changing group_size variable to alter the number of sequence files you want int each folder
#' If you have many cores to work with, you can split the sequences into many folders
#' The alignment_bash_script.txt file contains a simple script to align sequences. It assumes you have clustalo in your path, so make sure that is true
#' 
#' Complete the following actions on the command line:
#' 1) navigate to your fasta folder
#' 2) copy-paste the split_folders.txt contents into your command line environment and execute.
#' 3) start a screen instance, naming it after a specific folder (e.g. "screen -O dir001")
#' 4) navigate to a specific subdirectory (dir001, dir002, etc.) and execute the contents of the alignment_bash_script.txt
#' 5) exit the screen (screen + ctrl  + ctrl D) and repeat steps 3-4 for all subdirectories
#' After finishing the alignments, return to this R script


# calculate genetic summary stats
# calculate pi statistics for each species within each cell
aligned_folder <- "data/bold_seqs_aligned_100"

pi_df_one <- list.files(aligned_folder, full.names = TRUE) %>% 
  map_dfr(gen_dist_calc) 
###Write this to a csv. Can read in later if you don't want to run all of the stats again.
fwrite(pi_df_one, file = paste0(todays_results, "/output/pi_df_one_100.csv"))

print("Calculated pi and sequence statistics")

