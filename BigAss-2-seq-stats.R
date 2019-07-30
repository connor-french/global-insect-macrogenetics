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
library(tidyverse)
library(tidylog)
#library(furrr)
library(leaflet)
library(sf)
library(raster)
library(ape)
library(gdata)
library(entropy)

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
#' 
#' Calculate genetic diversity statistics for each cell. 
## ------------------------------------------------------------------------

#define output folder
#create directory to put fasta files
aligned_folder <- paste0(todays_results, "/output/bold_seqs_aligned")
dir.create(aligned_folder)

#split list into chunks for my janky "parallelization" (clustalo's threads option won't work on the lab computer)
#fasta_split <- split(list.files(fasta_folder), ceiling(seq_along(list.files(fasta_folder))/10))
fasta_files <- list.files(fasta_folder)

#align the fastas. Assumes that clustalo is in your path! if not, provide with the path in the argument
#choose your number of threads

map(fasta_files, 
    ~ system2("clustalo", 
            args =  c(paste0("-i ", getwd(), "/", fasta_folder, "/", .), paste0("-o ", getwd(), "/", aligned_folder, "/", .))
    )
)

#remove the folder with unaligned sequences
unlink(fasta_folder, recursive = TRUE)

#calculate genetic summary stats
#calculate pi statistics for each species within each cell
pi_df_one <- list.files(aligned_folder, full.names = TRUE) %>% 
  map_dfr(gen_dist_calc) 
###Write this to a csv. Can read in later if you don't want to run all of the stats again.
fwrite(pi_df_one, file = paste0(todays_results, "/output/pi_df_one.csv"))

print("Calculated pi and sequence statistics")

