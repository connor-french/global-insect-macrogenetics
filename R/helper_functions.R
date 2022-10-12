# helper functions for various scripts


##### %notin% ############
# a function to make it easier and more succinct to negate a vector for filtering
`%notin%` <- Negate(`%in%`)


###### hill_calc ########

## Get one hill number from a list of genetic distances. Original python code written by Isaac Overcast, with slight modifications (correct = TRUE implemented by CMF)
## dists are the OTU Tajima's pi
## order is the q order of the Hill number
## correct indicates if you want to correct for species richness or not. Default is TRUE
hill_calc <- function(dists, order = 1, correct = TRUE) { 
  if (order == 0) {
    return(length(dists))
  }
  if (order == 1) {
    h1 = exp(entropy::entropy(dists))
    if (correct) {
      return(h1 / length(dists))
    } else return(h1)
    
  }
  else {
    tot = sum(dists)
    proportions = dists/tot
    prop_order = proportions**order
    h2 = sum(prop_order)**(1/(1-order))
    if (correct) {
      return(h2 / length(dists))
    } else return(h2)
  }
}

###### Normalize predictors ########

normalize_vars <- function(df_in) {
  df_tib <- as_tibble(df_in) %>% 
    select(-geometry)
  
  df_sf <- df_in %>% 
    select(cell)
  
  norm_rec <- recipe(hill_1 ~ ., data = df_tib) %>%
    step_normalize(all_numeric(), -cell, -num_otu, -num_ind, -num_order, -contains("_pi"), -contains("hill_"))
  
  df_trans <- norm_rec %>%
    prep(training = df_tib) %>%
    juice(all_predictors()) %>%
    mutate(hill_1 = df_tib$hill_1)
  
  df_trans_sf <- left_join(df_sf, df_trans, by = "cell")
}

##### Join predictors with their respective predictors ######
join_predictors <- function(gen_sumstats, pred_rasts, resolution){
  exp_df <- pred_rasts[gen_sumstats$cell] %>% 
    as_tibble() %>% 
    mutate(cell = gen_sumstats$cell)
  
  if (resolution == "high") {
    full_df <- left_join(gen_sumstats, exp_df, by = "cell") %>% 
      mutate(soil_high_hwsd = as.factor(soil_high_hwsd),
             resolution = "high") %>% 
      arrange(cell)
  } else if (resolution == "medium") {
    full_df <- left_join(gen_sumstats, exp_df, by = "cell") %>% 
      mutate(soil_medium_hwsd = as.factor(soil_medium_hwsd),
             resolution = "medium") %>% 
      arrange(cell) 
  } else {
    full_df <- left_join(gen_sumstats, exp_df, by = "cell") %>% 
      mutate(resolution = "low") %>% 
      arrange(cell) 
  }
  
  
  return(full_df)
}

###### Filter dfs, removing plant phylo and missing data #####
filter_dfs <- function(df, resolution) {
  plant_res <- paste0("plant_", resolution, "_plant_phylo")
  df_filtered <- df %>% 
    select(-c(plant_res)) %>% 
    remove_missing()
  return(df_filtered)
}

##### sample_response_posterior ##########

# get draws from the posterior distribution
# response can be either gde or gdm, since those are what we're interested in
sample_response_posterior <- function(model, num_iterations = 1000, response = "gde") {
  
  # sample from the posterior distribution
  post_draws <- posterior_predict(model, iter = num_iterations) %>% 
    t() 
  colnames(post_draws) <- paste0("draw_", 1:num_iterations)
  
  d <- model$data
  
  # get the observed response vector
  if (response == "gde") {
    resp = d$gde
  } else if (response == "gdm") {
    resp = d$gdm
  }
  
  # create a data frame of posterior draws and add in the observed values for comparison.
  # in addition, I'm adding in the latitude and longitude 
  post_df <- post_draws %>% 
    as_tibble() %>% 
    mutate(observed = resp,
           longitude = model$data$lon_scaled * 1e6,
           latitude = model$data$lat_scaled * 1e6,
           id = 1:nrow(post_draws)) %>% 
    pivot_longer(cols = c(contains("_"), observed),
                 names_to = "draw",
                 values_to = "response_post") %>% 
    mutate(post_samples = if_else(draw == "observed", "observed", "posterior"))
  
  return(post_df)
}

##### tidy_post ##########

# function to retrieve the parameter posteriors and log probability posterior. This will also convert each class of posterior to a data frame with reasonable column names
tidy_post <- function(model) {
  post_dfs <- rstan::extract(model$model, permute = TRUE) %>% 
    map(as_tibble)
  
  # rename beta columns to their variable names
  colnames(post_dfs$B) <- colnames(model$X) %>% 
    janitor::make_clean_names(case = "snake")
  
  # rename spatial effect columns
  suffix <- str_remove(colnames(post_dfs$spatialEffectsKnots), "1.")
  colnames(post_dfs$spatialEffectsKnots) <- paste0("knot_", suffix)
  
  return(post_dfs)
}


######### resample_equal_area ############

# function to resample rasters to any equal area grid of the three resolutions we're considering.


# crs: define the crs you want to project to. It must be an equal area projection. The default crs is behrmann equal area. 
# x: a raster or raster stack
# km: km resolution you want (area will be km x km)
# is_categorical: if is_categorical is true, using "ngb" interpolation
resample_equal_area <- function(x, 
                                crs = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs", 
                                km,
                                is_categorical = FALSE) {
  
  # project raster to new coordinate system (this has no values- we put values back in later)
  x_ext <- raster::projectExtent(x, crs)
  
  # make the resolution square since we're interested in an equal area projection
  raster::res(x_ext) <- raster::xres(x_ext)
  
  # add original values to our projection. We can't take shortcuts like doing this instead of doing projectExtent first because it removes cells at the extreme latitudes.
  x_repro <- raster::projectRaster(x, x_ext)
  
  if (km == 96.5) {
    end_raster = template_high
  } else if (km == 193) {
    end_raster = template_medium
  } else if (km == 385.9) {
    end_raster = template_low
  } else stop("Incorrect resolution specification. Must be either 96.5 or 193")
  
  # resample our raster of interest to the beginning raster  
  if (is_categorical) {
    m <- "ngb"
  } else m <- "bilinear"
  
  
  x_resampled <- raster::resample(x_repro, end_raster, method = m) 
  
  
  return(x_resampled)
}


# Function from Gelman et al. 2018 to calculate a bayes version of the R2 value. I modified it to account for the Gaussian Process variance.  

bayes_R2_glmmfields <- function(fit) {
  y_pred <- glmmfields::posterior_linpred(fit)
  var_fit <- apply(y_pred, 1, var)
  var_res <- as.matrix(fit$model, pars = c("sigma"))^2
  var_gp <- as.matrix(fit$model, pars = c("gp_sigma"))^2
  return(var_fit / (var_fit + var_res + var_gp))
}

# Function to filter data sets for number of orders

order_filter <- function(order_name, analysis_data, min_n = 5) {
  raw_pi %>% 
    filter(order == order_name) %>% 
    group_by(cell) %>% 
    filter(n() >= 5, 
           cell %in% analysis_data$cell) %>% 
    summarize(gdm = sqrt(mean(pi)),
              gde = hill_calc(pi))
}


# function to run glmmfields models
run_glmmfields <- function(response, 
                           predictors, 
                           hypothesis, 
                           min_otu = 150, 
                           knots = 20, 
                           save_file = TRUE,
                           df, 
                           seed) {
 model <- glmmfields(
    formula(paste0(response," ~ ", predictors)),
    data = df,
    family = gaussian(),
    lat = "lat_scaled",
    lon = "lon_scaled",
    nknots = knots,
    iter = 7000,
    save_log_lik = TRUE,
    chains = 4,
    estimate_df = TRUE,
    covariance = "exponential",
    # intercept can't move beyond -1 or 1, so a relatively small scale is justified.
    prior_intercept = student_t(
      df = 1000,
      location = 0,
      scale = 1
    ),
    # betas are going to be very small too (definitely under 0.05), because the predictors are normalized and centered, and the response is bounded between zero and one (our data is between 0.35ish and 0.65ish, and Hill numbers wouldn't realistically reach the extremes). So a sigma of 0.1 is weakly regularizing and a normal prior is appropriate. 
    prior_beta = student_t(1000, 0, 0.1),
    prior_sigma = half_t(1000, 0, 1),
    prior_gp_theta = half_t(1000, 0, 5),
    prior_gp_sigma = half_t(1000, 0, 1),
    control = list(adapt_delta = 0.99,
                   max_treedepth = 15),
    seed = seed
  )
 
 hypo_file <- hypothesis %>% 
   paste0("-", response, "-", min_otu, ".rds")
 if (save_file) {
   write_rds(model, here("output", "models", hypo_file))
 }
 
 
 return(model)
 
}
  
  
# function to plot ELPD results for LOO model selection

plot_loo <- function(loo_output, sumstat, min_otu) {
  loo_output %>% 
    as_tibble(rownames = "hypothesis") %>% 
    mutate(hypothesis = fct_reorder(hypothesis, elpd_diff, .desc = TRUE)) %>% 
    ggplot() +
    geom_hline(yintercept = 0, color = "red") +
    geom_point(aes(x = hypothesis, y = elpd_diff)) +
    geom_errorbar(aes(x = hypothesis, y = elpd_diff,
                      ymin=elpd_diff - se_diff, ymax=elpd_diff + se_diff)) +
    labs(
      title = paste0(sumstat, ", Min OTU = ", min_otu),
      y = "Diff in ELPD") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20))
}

# function for order maps
map_order <- function(df, order_name, sumstat = "gdm", scale_guide = TRUE) {
  
  if (sumstat == "gdm") {
    # plot a red line for the median of the observed data
    obs_median <- median(full_sf_100$gdm)
    
    # for density plot limits and color scale limits
    lims <- c(min(bind_rows(outlier_df_list_sf)$gdm), max(bind_rows(outlier_df_list_sf)$gdm))
    
    # for plotting
    ss <- sym("gdm")
    
    fill_color <-  "mako"
    
  } else if (sumstat == "gde") {
    obs_median <- median(full_sf_100$gde)
    
    # for density plot limits and color scale limits
    lims <- c(min(bind_rows(outlier_df_list_sf)$gde), max(bind_rows(outlier_df_list_sf)$gde))
    
    # for plotting
    ss <- sym("gde")
    
    fill_color <- "rocket"
  }
  
  
  
  
  # num_otus for title
  if (order_name == "Observed data") {
    n_otus <- length(unique(pi_ordered$bin_uri))
  } else {
    n_otus <- length(unique(pi_outlier[pi_outlier$order == order_name,]$bin_uri))
  }
  
  
  # inset GDM histograms
  dens_inset <- ggplotGrob(ggplot(data = df, aes(x = !!ss)) +
                             geom_density(color = "black", fill = "gray") +
                             geom_density(data = full_sf_100, aes(x = !!ss), color = "black", fill = "transparent", linetype = "dashed") +
                             #geom_vline(xintercept = obs_median, color = "red") +
                             theme_insects() +
                             labs(x = str_to_upper(sumstat)) +
                             scale_x_continuous(limits = lims, labels = label_number(accuracy = 0.01)) +
                             theme(axis.title = element_text(size = 12),
                                   axis.text = element_text(size = 8),
                                   axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank()))
  
  if (scale_guide == TRUE) {
    ggplot() +
      geom_sf(data = world_base_map, fill = "lightgray", color = "lightgray") +
      geom_sf(data = df, aes(fill = !!ss, color = !!ss)) +
      # set the minimum and maximum value as the minimum and maximum across all data sets
      scale_fill_viridis_c(option = fill_color, 
                           direction = 1,
                           limits = lims) +
      scale_color_viridis_c(option = fill_color,
                            direction = 1,
                            limits = lims,
                            guide = "none") +
      labs(title = str_c(order_name, ", N=", nrow(df), 
                         " cells. N=", n_otus, " OTUs."),
           fill = str_to_upper(sumstat))  +
      theme_insects() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5,
                                      size = 25)) +
      annotation_custom(grob = dens_inset, 
                        xmin = -19000000, 
                        xmax = -10000000,
                        ymin = -7000000, 
                        ymax = 0)
  }
  else {
    ggplot() +
      geom_sf(data = world_base_map, fill = "lightgray", color = "lightgray") +
      geom_sf(data = df, aes(fill = !!ss, color = !!ss)) +
      scale_fill_viridis_c(option = fill_color, 
                           direction = 1) +
      scale_color_viridis_c(option = fill_color, 
                            direction = 1,
                            guide = "none") +
      labs(title = str_c(order_name, ", N=", nrow(df), " cells. N=", n_otus, " OTUs."),
           fill = str_to_upper(sumstat))  +
      theme_insects() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, 
                                      size = 25)) +
      annotation_custom(grob = dens_inset, 
                        xmin = -19000000, 
                        xmax = -10000000,
                        ymin = -7000000, 
                        ymax = 0)
  }
  
}


# tidy t-test output
tidy_ttest <- function(x, y, order, sumstat) {
  t <- t.test(x, y)
  mean_x <- t$estimate["mean of x"]
  mean_y <- t$estimate["mean of y"]
  upper_95 <- t$conf.int[1]
  lower_95 <- t$conf.int[2]
  t_stat <- t$statistic
  df <- t$parameter
  p <- round(t$p.value, 3)
  
  out_tib <- tibble(
    mean_x = mean_x,
    mean_y = mean_y,
    diff_in_means = mean_x - mean_y,
    upper_95 = upper_95,
    lower_95 = lower_95,
    t_statistic = t_stat,
    df = df,
    p_value = p,
    order = order,
    sumstat = sumstat
  )
  
  return(out_tib)
}


# tidy correlation output
tidy_cor <- function(x, y, var, sumstat) {
  t <- cor.test(x, y)
  rho <- t$estimate
  upper_95 <- t$conf.int[1]
  lower_95 <- t$conf.int[2]
  t_stat <- t$statistic
  df <- t$parameter
  p <- round(t$p.value, 3)
  
  out_tib <- tibble(
    rho = rho,
    upper_95 = upper_95,
    lower_95 = lower_95,
    t_statistic = t_stat,
    df = df,
    p_value = p,
    var = var,
    sumstat = sumstat
  )
  
  return(out_tib)
}

# plot partial dependence curves

plot_part_dep <- function(df, variable) {
  
  if (variable == "current_medium_bio_13") {
    x_lab <- "WET"
  } else if (variable == "current_medium_bio_15") {
    x_lab <- "Precip. seasonality"
  } else if (variable == "current_medium_bio_2") {
    x_lab <- "Temp. seasonality"
  } else if (variable == "current_medium_bio_5") {
    x_lab <- "WARM"
  } else stop("Choose a variable from the final model.")
  
  df_name <- deparse(substitute(df))
  
  if (df_name == "beta_posts_gde") {
    y_lab <- "GDE"
    response <- sym("gde")
  } else if (df_name == "beta_posts_gdm") {
    y_lab <- "GDM"
    response <- sym("gdm")
  } else stop("The input data frame needs to be either beta_posts_gde or beta_posts_gdm.")
  
  var_name <- sym(variable)
  
  
  
  ggplot(data = model_data, aes(x = {{ var_name }}, y = {{ response }})) +
    # to easily set up the plot dimensions, making scatterplot without the points
    geom_point(color = "transparent") +
    geom_abline(data = df, aes(intercept = intercept, slope = {{ var_name }}), alpha = 0.2, color = "darkgray") +
    geom_abline(intercept = median(df$intercept), slope = median(df[[var_name]]), color = "darkgreen", size = 1.25) +
    labs(x = x_lab,
         y = y_lab) +
    theme_insects()
  
}


# function to help with predictor maps
plot_predictor <- function(variable) {
  pred_plot <- 
    ggplot() +
    geom_sf(data = all_predictors, aes(fill = .data[[variable]], color = .data[[variable]])) +
    # scale_fill_gradientn(colors = pal) +
    # scale_color_gradientn(colors = pal, guide = NULL) + 
    scale_fill_viridis_c() +
    scale_color_viridis_c(guide = NULL) +
    geom_sf(data = world_base_coast, fill = "transparent") +
    labs(fill = variable) +
    theme_insects() +
    theme(legend.text = element_text(size = 10),
          legend.title = element_text(size = 14))
  
  return(pred_plot)
  
}
