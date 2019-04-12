sumstat_plot_fun <- function(data, stat){
  ggplot(data, aes_string(x = stat, y = "..density..")) +
    geom_density(fill = "grey", alpha = 0.4) +
    geom_histogram(fill = "darkorange", color = "black", alpha = 0.7) +
    labs(y = "Density") +
    theme_minimal() 
}