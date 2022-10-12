# perform pca on the predictors to make sure the predictor space in the validation data isn't outside the range of the training data 
# extrapolation would be indicated by test data being clustered at the edges of the PCA plot
run_predictor_pca <- function(train_data, folds) {
  recipe(~., data = train_data) %>% 
    step_pca(all_predictors()) %>% 
    prep() %>% 
    bake(new_data = NULL) %>%
    bind_cols(train_data) %>% 
    mutate(split_1 = if_else(row_number() %in% folds$splits[[1]]$in_id, "train", "test"),
           split_2 = if_else(row_number() %in% folds$splits[[2]]$in_id, "train", "test"),
           split_3 = if_else(row_number() %in% folds$splits[[3]]$in_id, "train", "test"),
           split_4 = if_else(row_number() %in% folds$splits[[4]]$in_id, "train", "test"),
           split_5 = if_else(row_number() %in% folds$splits[[5]]$in_id, "train", "test")) %>% 
    pivot_longer(starts_with("split"), names_to = "split_number", values_to = "cv_split")
}



# convert fold data to tibbles
# for some reason, tidymodels can't handle the sf data in the folds

ssf_to_sf <- function(splits) {
  for (i in 1:length(splits)) {
    splits[[i]]$data <- as_tibble(splits[[i]]$data) %>% select(-geometry)
  }
  return(splits)
}


# center and scale data
scale_data <- function(df) {
  c <- st_coordinates(st_centroid(df))
  df_ll <- df %>%
    as_tibble() %>%
    select(-geometry) %>%
    mutate(
      lon = c[, "X"],
      lat = c[, "Y"],
      lon_scaled = 1e-5 * lon,
      lat_scaled = 1e-5 * lat
    )
  
  pred_vars <- c("bio_2", 
                 "bio_5", 
                 "bio_13", 
                 "bio_14", 
                 "bio_15", 
                 "ghh_std_dev", 
                 "human_gHM", 
                 "precip_trend", 
                 "precip_var",
                 "temp_trend",
                 "temp_var")
  
  df_scaled <- df_ll %>% 
    mutate(across(all_of(pred_vars), ~scale(.x)[,1]))
  
  return(df_scaled)
}


run_stan <- function(response = "gde", df) {
  
  pred_vars <- c("bio_2", 
                 "bio_5", 
                 "bio_13", 
                 "bio_14", 
                 "bio_15", 
                 "ghh_std_dev", 
                 "human_gHM", 
                 "precip_trend", 
                 "precip_var",
                 "temp_trend",
                 "temp_var")
  
  pred_vars <- paste(pred_vars, collapse = " + ")
  
  m_form <- formula(
    paste(
      response,
      "~",
      pred_vars
      )
  )
  
  mod <- stan_glm(
    m_form,
    data = df,
    prior = normal(location = 0, scale = 0.1),
    prior_intercept = normal(location = 0, scale = 1),
    family = gaussian(),
    iter = 4000,
    chains = 4
  )
  
  return(mod)
  
}


# projpred selection
run_projpred <- function(mod) {
  ref <- get_refmodel(mod)
  
  cv <- cv_varsel(ref,
                  method = "forward")
  return(cv)
}


# model glmmfields
run_glmmfields <- function(response = "gde", df, pp, nknots = 20) {
  
  # get number of vars to keep
  var_size <- projpred::suggest_size(pp)
  
  # get names of vars
  var_list <- projpred::solution_terms(pp)
  
  # subset names for "keeper" variables
  var_keep <- var_list[1:var_size]
  
  # format for modeling
  pred_vars <- paste(var_keep, collapse = " + ")
  
  m_form <- formula(
    paste(
      response,
      "~",
      pred_vars
    )
  )
  
  glmmfields(
    m_form,
    data = df,
    family = gaussian(),
    lat = "lat_scaled",
    lon = "lon_scaled",
    nknots = nknots,
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
                   max_treedepth = 15)
  )
} 


# function to get the posteriors for the model coefficients of each model
get_beta_posts <- function(mod_obj) {
  # only get post-warmup draws
  posts <- as.matrix(mod_obj$model)[3501:7000,]
  # extract beta coefficients
  posts <- posts[, grep("^B", colnames(posts))]
  
  if (ncol(posts) > 2) {
    # remove intercept
    posts <- posts[,-1]
    
    # get real variable names
    colnames(posts) <- colnames(mod_obj$X)[-1]
  } else {
    colnames(posts) <- colnames(mod_obj$X)
    
    # assigning the intercept dummy values so I can see the coefficient posteriors
    posts[,"(Intercept)"] <- posts[,2]
    
    colnames(posts)[1] <- "Dummy"
    
  }
  
  return(posts)
}

# I have to scale the test data by the same mean and sd I scaled the training data
scale_test <- function(df_train, df_test) {
  c <- st_coordinates(st_centroid(df_test))
  
  df_test <- df_test %>%
    as_tibble() %>%
    select(-geometry) %>%
    mutate(
      lon = c[, "X"],
      lat = c[, "Y"],
      lon_scaled = 1e-5 * lon,
      lat_scaled = 1e-5 * lat
    )
  
  pred_vars <- c("bio_2", 
                 "bio_5", 
                 "bio_13", 
                 "bio_14", 
                 "bio_15", 
                 "ghh_std_dev", 
                 "human_gHM", 
                 "precip_trend", 
                 "precip_var",
                 "temp_trend",
                 "temp_var")
  
  for (i in pred_vars) {
    mean_train <- mean(df_train[[i]])
    sd_train <- sd(df_train[[i]])
    df_test[, i] <- (df_test[[i]] - mean_train) / sd_train
  }
  
  
  
  return(df_test)
}


# observed vs predicted plotting
plot_obs_vs_pred <- function(df, resp) {
  
  lt <- summary(lm(formula(paste(resp, "~", ".pred")), data = df))
  
  r_2 <- lt$r.squared %>% 
    round(2)
  
  slope <- lt$coefficients[,1][[".pred"]] %>% 
    round(2)
  
  intercept <- lt$coefficients[,1][["(Intercept)"]] %>% 
    round(3)
  
  rmse <- yardstick::rmse_vec(truth = df[[resp]], estimate = df[[".pred"]]) %>% 
    round(3)
  
  ggplot(df, aes_string(x = ".pred", y = resp)) +
    geom_point(fill = "darkgray", color = "black", shape = 21) +
    geom_smooth(color = "red",
                method = "lm",
                formula = "y ~ x") +
    geom_abline(slope = 1, intercept = 0) +
    labs(caption = str_wrap(paste0("RMSE = ", rmse, ", R^2 = ", r_2, ", Slope = ", slope, ", Y-int = ", intercept)), 
         width = 10,
         x = "Predicted", y = "Observed") +
    theme_bw()
}

# function to map test predictions and error

map_test <- function(df, var = c(".pred", ".pred_ci", ".resid", ".log_ci")) {
  ggplot() +
    geom_sf(data = st_geometry(world_small)) +
    geom_sf(data = df, aes_string(color = var[1], fill = var[1])) +
    scale_color_viridis_c() +
    scale_fill_viridis_c() +
    theme_bw()
}



run_latitude <- function(response = "gde", df, nknots = 20) {
  

  m_form <- formula(
    paste(
      response,
      "~ abs(lat_scaled)"
    )
  )
  
  glmmfields(
    m_form,
    data = df,
    family = gaussian(),
    lat = "lat_scaled",
    lon = "lon_scaled",
    nknots = nknots,
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
                   max_treedepth = 15)
  )
} 

  

# function to rasterize sf polygons of predictor variables for MESS analysis

rasterize_preds <- function(df, rast_template, pred_vars) {
  rasters <- map(pred_vars, ~fasterize::fasterize(sf = df, raster = rast_template, field = .x))
  s <- stack(rasters) %>% round(5)
  names(s) <- pred_vars
  return(s)
}

# function to convert the mess values to a polygon layer for easy masking of the final maps
mess_to_poly <- function(mess) {
  mess %>%
    rasterToPolygons() %>%
    st_as_sf(crs = crs_behr) %>%
    mutate(
      mess_val = na_if(mess, "Inf"),
      mess_binary = if_else(mess_val < -1e-8, "non_analogous", "analogous")
    ) %>%
    select(-mess)
}

# make maps of mess models

map_mess <- function(mess) {
  ggplot() +
    geom_sf(data = st_geometry(world_small)) +
    geom_sf(data = mess, aes(color = mess_binary, fill = mess_binary)) +
    scale_color_manual(values = c("#56B4E9", "#E69F00"), na.value = "transparent", aesthetics = c("color", "fill"), na.translate = FALSE) +
    theme_bw()
}


