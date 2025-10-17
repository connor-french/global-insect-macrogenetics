library(tidyverse)
library(here)
library(patchwork)

bold_cell_sum <- read_csv(here("output", "spreadsheets", "bold_cell_summaries.csv")) %>% 
  mutate(otu_filter = as.factor(otu_filter),
         res_filter = fct_relevel(res_filter, "high", "medium", "low")) 

bold_filt_sum <- read_csv(here("output", "spreadsheets", "bold_summary_resolutions.csv")) %>% 
  mutate(otu_filter = as.factor(otu_filter),
         res_filter = fct_relevel(res_filter, "high", "medium", "low")) 

bold_filt_allsum <- left_join(bold_cell_sum, 
                              bold_filt_sum, 
                              by = c("res_filter", "ind_filter", "otu_filter"))


bold_filt_sum %>% 
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






