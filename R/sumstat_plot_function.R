sumstat_plot_fun <- function(data, stat){
  ggplot(data, aes_string(x = stat)) +
    geom_histogram(fill = "darkorange", color = "black", alpha = 0.7, bins = 50) +
    labs(y = "Frequency") +
    theme_minimal() 
}