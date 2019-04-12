#function to calculate pi for each cell number. I'm using dplyr::summarise because plyr also has a summarise function that can create a lot of confusion. nspec is the number of species per cell that were filtered
pi_summary_fun <- function(df, n_species) {
  gdf <- group_by(df, cells)
  sdf <- sample_n(gdf, n_species, replace = TRUE)
  sum_df <- dplyr::summarise(sdf, median.pi = median(avg_pi), mean.pi = mean(avg_pi), sd.pi = sd(avg_pi), hill.zero = hill_calc(avg_pi, 0), hill.one = hill_calc(avg_pi, 1), hill.two = hill_calc(avg_pi, 2), shannon = entropy::entropy(avg_pi))
  sum_df
}