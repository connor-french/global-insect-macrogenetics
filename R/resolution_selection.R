library(tidyverse)
library(here)
library(patchwork)

pw_pi_low <- read_csv(here("output", "spreadsheets", "cell_low_3_10_pi.csv")) 
pw_pi_med <- read_csv(here("output", "spreadsheets", "cell_medium_3_10_pi.csv"))
pw_pi_hi <- read_csv(here("output", "spreadsheets", "cell_high_3_10_pi.csv"))


res <- c(min_10 = 10, min_25 = 25, min_50 = 50, min_100 = 100, min_150 = 150, min_200 = 200)

summarize_sampling <- function(pi_df, res) {
  pdf <- pi_df %>% 
    group_by(cell) %>% 
    summarize(num_otu = length(bin_uri),
              num_ind = sum(num_ind),
              num_order = length(unique(order)))
  
  pi_sum <- map_df(res, ~filter(pdf, num_otu >= .x), .id = "min_otu")
  
  return(pi_sum)
}

pi_sum_low <- summarize_sampling(pi_df = pw_pi_low, res = res) %>% 
  mutate(resolution = "low")
pi_sum_med <- summarize_sampling(pi_df = pw_pi_med, res = res) %>% 
  mutate(resolution = "med")
pi_sum_hi <- summarize_sampling(pi_df = pw_pi_hi, res = res) %>% 
  mutate(resolution = "hi")


pi_sum_all <- bind_rows(pi_sum_low, pi_sum_med, pi_sum_hi)

cell_sum_all <- pi_sum_all %>% 
  group_by(resolution, min_otu) %>% 
  summarize(num_cells = length(cell),
            median_otu = median(num_otu),
            total_otu = sum(num_otu),
            median_ind = median(num_ind),
            total_ind = sum(num_otu),
            median_order = median(num_order),
            total_order = sum(num_order),
            var_otu = var(num_otu)) %>% 
  ungroup() %>% 
  mutate(min_otu_num = as.numeric(str_remove_all(min_otu, "min_")),
         min_otu = fct_reorder(min_otu, min_otu_num),
         resolution = fct_relevel(resolution, levels = c("hi", "med", "low")))

num_cells_plot <- cell_sum_all %>% 
  ggplot(aes(x = resolution, y = num_cells, fill = min_otu), color = "black") + 
  geom_point(shape = 21, size = 3)  +
  scale_fill_viridis_d()

median_order_plot <- cell_sum_all %>% 
  ggplot(aes(x = resolution, y = median_order, fill = min_otu), color = "black") + 
  geom_point(shape = 21, size = 3)  +
  scale_fill_viridis_d()

median_otu_plot <- cell_sum_all  %>% 
  ggplot(aes(x = resolution, y = median_otu, fill = min_otu), color = "black") + 
  geom_point(shape = 21, size = 3)  +
  scale_fill_viridis_d()

median_ind_plot <- cell_sum_all %>% 
  ggplot(aes(x = resolution, y = median_ind, fill = min_otu), color = "black") + 
  geom_point(shape = 21, size = 3)  +
  scale_fill_viridis_d()

total_ind_plot <- cell_sum_all %>% 
  ggplot(aes(x = resolution, y = total_ind, fill = min_otu), color = "black") + 
  geom_point(shape = 21, size = 3)  +
  scale_fill_viridis_d()

total_otu_plot <- cell_sum_all %>% 
  ggplot(aes(x = resolution, y = total_otu, fill = min_otu), color = "black") + 
  geom_point(shape = 21, size = 3)  +
  scale_fill_viridis_d()

var_otu_plot <- cell_sum_all %>% 
  ggplot(aes(x = resolution, y = var_otu, fill = min_otu), color = "black") + 
  geom_point(shape = 21, size = 3)  +
  scale_fill_viridis_d()


total_res_plot <- num_cells_plot + median_order_plot +
  var_otu_plot + median_otu_plot + total_otu_plot + 
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")
  
total_res_plot

ggsave(filename = here("output", "publication_figs", "resolution_decision_fig.png"),
       plot = total_res_plot,
       width = 36,
       height = 22.5,
       dpi = 300,
       unit = "cm")


epi_sum_all %>% 
  mutate(ind_filter = as.factor(ind_filter),
         otu_filter = as.factor(otu_filter))  %>% 
  ggplot(aes(x = otu_median, 
             y = otu_var, 
             color = num_cells, 
             shape = res_filter)) + 
  geom_point(size = 5) +
  geom_text(aes(label = paste0("Ind: ", ind_filter, ", ", "OTU: ", otu_filter)), nudge_x = 0.5, nudge_y = -1500) +
  geom_vline(xintercept = median(bold_filt_sum$otu_median)) +
  geom_hline(yintercept = median(bold_filt_sum$otu_var)) +
  scale_color_distiller(palette = "PRGn") +
  #xlim(4, 12) +
  theme_dark()



# Number of orders per cell for all filter regimes
# Looks like the min_otu = 100 is consistently low for all resolution levels
bold_cell_sum %>% 
  filter(ind_filter == 3) %>% 
  ggplot(aes(x = n_order, color = otu_filter)) +
  geom_boxplot() +
  facet_wrap(~res_filter, nrow = 3) +
  theme_minimal()

# number of OTUs per cell for each filter
bold_cell_sum %>% 
  filter(ind_filter == 3) %>% 
  ggplot(aes(x = n_otu, y = otu_filter, color = )) +
  geom_boxplot() +
  xlim(NA, 4000) +
  facet_wrap(~res_filter, nrow = 3) +
  theme_minimal()

  

bold_all_sum <- bold_cell_sum %>% 
  filter(ind_filter == 3) %>% 
  group_by(res_filter, ind_filter, otu_filter) %>% 
  summarize(total_otu = sum(n_otu),
            median_ind = median(total_ind),
            total_ind = sum(total_ind),
            median_otu = median(n_otu),
            median_order = median(n_order),
            n_cells = n(),
            lower_90 = HDInterval::hdi(n_otu, credMass = 0.9)[1],
            upper_90 = HDInterval::hdi(n_otu, credMass = 0.9)[2]) %>% 
  mutate(dif_hdi = upper_90 - lower_90) %>% 
  ungroup()

# variation in num_OTUs per cell across cells ~ median number of OTUs per cell
median_order <- bold_all_sum %>% 
  filter(ind_filter == 3) %>% 
  mutate(pt_label = str_c(res_filter, otu_filter, sep = "_")) %>% 
  ggplot(aes(x = median_otu, y = n_cells, label = pt_label)) +
  geom_point() +
  geom_label(aes(color = median_order),
             size = 2) + 
  labs(x = NULL, y = NULL,
    color = "Median # orders") +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_color_viridis_c() + 
  theme_minimal()

variation_otus <- bold_all_sum %>% 
  filter(ind_filter == 3) %>% 
  mutate(pt_label = str_c(res_filter, otu_filter, sep = "_")) %>% 
  ggplot(aes(x = median_otu, y = n_cells, label = pt_label)) +
  geom_label(aes(color = dif_hdi),
             size = 2) + 
  labs(x = NULL, y = NULL,
    color = "Variation in # OTUs") +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_color_viridis_c() +
  theme_minimal()

# total number of individuals
total_inds <- bold_all_sum %>% 
  filter(ind_filter == 3) %>% 
  mutate(pt_label = str_c(res_filter, otu_filter, sep = "_")) %>% 
  ggplot(aes(x = median_otu, y = n_cells, label = pt_label)) +
  geom_point() +
  geom_label(aes(color = total_ind),
             size = 2) + 
  labs(x = NULL, y = NULL,
    color = "Total # individuals") +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_color_viridis_c() +
  theme_minimal()

# total number of OTUs (not unique)
total_otus <- bold_all_sum %>% 
  filter(ind_filter == 3) %>% 
  mutate(pt_label = str_c(res_filter, otu_filter, sep = "_")) %>% 
  ggplot(aes(x = median_otu, y = n_cells, label = pt_label)) +
  geom_point() +
  geom_label(aes(color = total_otu), 
             size = 2) + 
  labs(x = NULL, y = NULL,
    color = "Total # OTUs") +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_color_viridis_c() +
  theme_minimal()

otu_meta_plot <- variation_otus + 
  median_order +
  total_otus + 
  plot_annotation(title = str_wrap("Number of cells as a function of median number OTUs per cell"),
                  tag_levels = "a",
                  tag_suffix = ")")+
  plot_layout(ncol = 2)

otu_meta_plot


ggsave(filename = here("output", "publication_figs", "resolution_decision_plot.svg"),
       plot = otu_meta_plot,
       width = 24,
       height = 15,
       dpi = 300,
       unit = "cm")






