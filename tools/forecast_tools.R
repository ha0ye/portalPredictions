"%>%" <- magrittr::"%>%"
"quos" <- rlang::"quos"
"!!!" <- rlang::"!!!"

#' Combine all new forecasts and AIC files (from the tmp directory) and add 
#'   ensembles
#' 
#' @param forecast_date
#' @param filename_suffix
#' @return list(forecasts, all_model_aic)
#'
forecastall <- function(forecast_date, filename_suffix = "forecasts"){
  
  forecast_files <- paste(filename_suffix, ".csv", sep = "")
  forecast_filelist <- list.files("tmp", forecast_files, full.names = TRUE)
  col_classes <- c("Date", "integer", "integer", "integer", "character", 
                   "character", "character", "character", "numeric", 
                   "numeric", "numeric", "integer", "integer", "integer")
  forecasts <- do.call(rbind, 
                 lapply(forecast_filelist, read.csv, na.string = "", 
                   colClasses = col_classes))

  aic_files <- paste(filename_suffix, "_model_aic.csv", sep = "")
  aic_filelist <- list.files("tmp", aic_files, full.names = TRUE)
  all_model_aic <- do.call(rbind,
                     lapply(aic_filelist, read.csv, na.strings = ""))

  fcast_date <- as.character(forecast_date)
  forecast_suffix <- paste(fcast_date, filename_suffix, ".csv", sep = "")
  forecast_filename <- file.path("predictions", forecast_suffix)
  aic_suffix <- paste(fcast_date, filename_suffix, "_model_aic.csv", sep = "")
  model_aic_filename <- file.path("predictions", aic_suffix)
  append_csv(forecasts, forecast_filename)
  append_csv(all_model_aic, model_aic_filename)
  
  ensemble <- make_ensemble(forecasts) %>% 
                subset(select = colnames(forecasts))
  append_csv(ensemble, forecast_filename)
  
  return(list(forecasts, all_model_aic))
}

#' Append a csv without re-writing the header
#'
#' @param df data frame
#' @param filename existing file
#'
append_csv <- function(df, filename){
  write.table(df, filename, sep = ',', append = file.exists(filename),
    row.names = FALSE, col.names = !file.exists(filename))
}

#' Get all model aic values and calculate akaike weights
#'
#' @param forecast_folder folder with the forecast files
#'
compile_aic_weights <- function(forecast_folder = "./predictions"){
  f_names <- list.files(forecast_folder, full.names = TRUE, recursive = TRUE)
  which_aic_files <- which(grepl('model_aic', f_names))
  aic_filenames <- f_names[which_aic_files]

  read_formula <- ~read.csv(.x, na.strings = '', stringsAsFactors = FALSE)
  all_aic <- purrr::map(aic_filenames, read_formula) %>% 
               dplyr::bind_rows()

  wts <- all_aic %>%
           dplyr::group_by(date, currency, level, species, fit_start_newmoon, 
             fit_end_newmoon, initial_newmoon) %>%
           dplyr::mutate(delta_aic = aic - min(aic), 
             weight = exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))) %>%
           dplyr::ungroup()
  return(wts)
}

#'Create the ensemble model from all other forecasts using the weighted
#'  mean and weighted sample variance
#'
#' @details   Mean is the weighted mean of all model means, variance is the
#'  weighted mean of all model variances + the variances of the weighted mean 
#'  using the unbiased estimate of sample variance. See 
#'  https://github.com/weecology/portalPredictions/pull/65 we only store the
#'  prediction interval for models, so backcalculate individual model 
#'  variance assuming the same CI_level throughout. 
#'
#' @param all_forecasts forecasts
#' @param models_to_use specific models to use
#' @param CI_level confidence interval level
#'
make_ensemble <- function(all_forecasts, models_to_use = NA, CI_level = 0.9){

  weights <- compile_aic_weights()
  weights$date <- as.Date(weights$date)
  CI_quantile <- qnorm((1 - CI_level) / 2, lower.tail = FALSE)

  lft_jn_cols <- c("date", "model", "currency", "level", "species", 
                   "fit_start_newmoon", "fit_end_newmoon", "initial_newmoon")
  
  model_var <- quos(model_var = ((UpperPI - estimate) / CI_quantile)^2)
  grouping <- quos(date, newmoonnumber, forecastmonth, forecastyear, level, 
                currency, species, fit_start_newmoon, fit_end_newmoon, 
                initial_newmoon)
  summarising <- quos(ensemble_estimate = sum(estimate * weight), 
                   weighted_ss = sum(weight * (estimate - ensemble_estimate)^2),
                   ensemble_var = sum(model_var * weight) + 
                                  weighted_ss / (n()*sum(weight)-1),
                   sum_weight = sum(weight))
  wtd_ests <- all_forecasts %>%
              dplyr::mutate(!!!model_var) %>%
              dplyr::left_join(weights, by = lft_jn_cols) %>%
              dplyr::group_by(!!!grouping) %>%
              dplyr::summarise(!!!summarising) %>% 
              dplyr::ungroup() 
              
  # Check that the summed weight of all the model ensembles is 1
  summed_weights <- round(wtd_ests$sum_weight, 10)
  if (!all(summed_weights == 1 | is.na(wtd_ests$sum_weight))){ 
    stop("Summed weights do not equal 1")
  }

  PI <- quos(LowerPI = ensemble_estimate - (sqrt(ensemble_var) * CI_quantile),
          UpperPI = ensemble_estimate + (sqrt(ensemble_var) * CI_quantile))
  ensemble <- wtd_ests %>%
              dplyr::mutate(!!!PI) %>%
              dplyr::mutate(LowerPI = ifelse(LowerPI<0, 0, LowerPI)) %>%
              dplyr::rename(estimate = ensemble_estimate) %>%
              dplyr::select(-ensemble_var, -weighted_ss, -sum_weight)
  
  ensemble$model <- "Ensemble"
  return(ensemble)
}


#' get species-level predictions
#'
#' @param data
#' @param lvl
#' @param lead_time
#'
get_sp_predicts <- function(data, lvl, lead_time){

  data <- transform(data, forecast_date = format(date, "%b %Y"))
  data <- transform(date = as.Date(date, "%Y-%m-%d"))
  data1 <- dplyr::filter(data, level == lvl, date == max(as.Date(date)))
  target_moon <- min(data1$newmoonnumber) + (lead_time - 1)
  data2 <- dplyr::filter(data1, newmoonnumber == target_moon)
}

#' this is the second plot on the "Species-level Forecast" page on the 
#'   "'Current Forecast" page on the website
#' 
#' 
#' @param data
#' @param title main title for plot
#' @return sp_predict is a plot object -- plot(sp_predict) displays it
#' 
plot_species_forecast <- function(data, title) {

  main_url <- "https://raw.githubusercontent.com/weecology/PortalData/master/"
  moon_suffix <- "Rodents/moon_dates.csv"
  moon_url <- paste(main_url, moon_suffix, sep = "")
  spp_suffix <- "Rodents/Portal_rodent_species.csv"
  spp_url <- paste(main_url, spp_suffix, sep = "")

  newmoons_table <- read.csv(text = RCurl::getURL(moon_url))
  target_nm <- unique(data$newmoonnumber)

  period_code <- newmoons_table %>%
                 dplyr::filter(newmoons_table$newmoonnumber == target_nm) %>%
                 dplyr::select(period) %>%
                 as.integer()

  spp_tbl_url <- paste("https://raw.githubusercontent.com/weecology/",
                   "PortalData/master/Rodents/Portal_rodent_species.csv", 
                   sep = "")
  species_table <- read.csv(text = RCurl::getURL(spp_url), 
                     stringsAsFactors = F, na.strings = "")
  species_names <- species_table %>% 
                   dplyr::select("speciescode", "scientificname") %>% 
                   rbind(c("total", "total")) %>%
                   merge(data[ , c("species", "estimate")], 
                     by.x = "speciescode", by.y = "species")
  
  sp_predict <- ggplot2::ggplot(data,
                  aes(x = estimate, y = reorder(species, estimate), 
                    xmin = LowerPI, xmax = UpperPI)) +
                ggplot2::geom_point() +
                ggplot2::geom_errorbarh() +
                ggplot2::ggtitle(title) + 
                ggplot2::ylab("Species") +
                ggplot2::xlab("Abundance") +
                ggplot2::scale_y_discrete(
                  breaks = reorder(data$species, data$estimate), 
                  labels = reorder(species_names$scientificname, 
                             species_names$estimate))
  
  return(sp_predict)
}

#' Compares forecasts to observations over different lead times.
#'
#' @details Error can be any function. The level, species, and currency 
#'   columns from observations and forecasts must have matching values. Note 
#'   this gives an average error value over many forecast iterations.
#'   Will only return values where there are matching comparison columns 
#'   (currency, level, species)
#'
#' @param observations dataframe Has the columns newmoonnumber, currency, 
#'   level, species, actual
#' @param forecasts dataframe passes the forecast validity check. Must have 
#'   matching values in the comparison columns
#' @param error_metric chr either 'RMSE' for root mean squared error or 
#'  "coverage" for the coverage of the prediction intervals
#' @param ci_value int The value of the forecast confidence interval to scale 
#'   PI values for the likelihood metric
#'
#' @return data.frame Data.frame with the columns model, error, lead_time, 
#'   level, species, currency
#'
calculate_forecast_error <- function(observations, forecasts, 
                                     error_metric = "RMSE", CI_level = 0.9){

  observations <- as.data.frame(observations)
  forecasts <- as.data.frame(forecasts)

  if (!forecast_is_valid(forecasts)){
    stop('Forecast dataframe not valid')
  }

  valid_obs_cols  <- c("newmoonnumber", "currency", "level", "species", 
                       "actual")
  if(!all(valid_obs_cols %in% colnames(observations))){
    stop("observation data.frame does not have valid column names")
  }

  # At least 1 matching value must be in each of these columns in the 
  # observations and forecasts
  # TODO: Ensure matching rows in all 3 columns at once instead of just one
  # at a time.
  column_check <- c()
  for(column in c("currency", "level", "species")){
    if(!any(unique(observations[ , column]) %in% unique(forecasts[ , column]))){
      column_check <- c(column_check, column)
    }
  }
  if(length(column_check) > 0){
    stop_msg <- paste("Columns do not match: ", column_check, collapse = " ")
    stop(stop_msg)
  }

  # Summarize to mean error by lead time. Lead time is number of new moons 
  # ahead of when the forecast was made. This assumes a forecast was made
  # with only the data available prior to the first NewMoonDate in the series.
  # TODO: Make the lead time the actual days or weeks once more frequent 
  #   forecasts are being made(see #37)
  forecasts <- forecasts %>%
                 dplyr::mutate(lead_time = newmoonnumber - initial_newmoon)
  
  # Calculate error
  in_jn_cols <- c("newmoonnumber", "currency", "level", "species")
  if (error_metric == "RMSE"){
    grouping <- quos(model, currency, level, species, lead_time)
    comparisons <- forecasts %>%
                   dplyr::inner_join(observations, by = in_jn_cols) %>%
                   dplyr::mutate(error_value = (estimate - actual)^2) %>%
                   dplyr::group_by(!!!grouping) %>%
                   dplyr::summarize(error_value = sqrt(mean(error_value))) %>%
                   dplyr::ungroup() %>%
                   dplyr::mutate(error_metric = "RMSE")

  } else if (error_metric == "coverage"){
    grouping <- quos(model, currency, level, species, lead_time, error_metric)
    new_vars <- quos(within_prediction_interval = actual >= LowerPI & 
                     actual <= UpperPI, error_metric = "coverage")
    summarising <- quos(error_value = mean(within_prediction_interval))
    comparisons <- forecasts %>%
                   dplyr::inner_join(observations, by = in_jn_cols) %>%
                   dplyr::mutate(!!!new_vars) %>%
                   dplyr::group_by(!!!grouping) %>%
                   dplyr::summarize(!!!summarising) %>%
                   dplyr::ungroup() %>%
                   dplyr::mutate(error_metric = "coverage")
  } else if (error_metric == "deviance"){
    stop("Deviance not implimented  yet")
  } else{
    stop(paste0("Error metric unknown: ", error_metric))
  }

  return(comparisons)
}

#' Plot the output of calculate_forecast_error(). Lead time on the x-axis,
#' error on the y-axis, different colored lines are different models.
#'
#' @param error_df data.frame The output from calculate_foreast_error()
#' @param level str Valid level
#' @param species str Valid species
#' @param currency str Valid currency
#' @param error_metric str error metric used
#'
plot_lead_time_errors <- function(error_df, level, species, currency, 
                           error_metric){
  plot_title <- paste0("Level: ", level, ", Species: ", species, ", 
                  Currency: ", currency)

  graph <- ggplot2::ggplot(error_df, aes(x = lead_time, y = error, 
             group = model, color = model)) +
           ggplot2::geom_point() +
           ggplot2::geom_line() +
           ggplot2::labs(y = error_metric, x = "Lead Time (New Moons)", 
             title = plot_title)
  plot(graph)
}

#' Ensure that a forecast file is in the correct format
#'
#' Tools for working with forecast data expect a certain format.
#' This ensures a forecast file meets those formats. All column
#' and variable names are case sensitive. For specification see:
#' https://github.com/weecology/portalPredictions/wiki/forecast-file-format
#'
#' @param forecast_df dataframe A dataframe read from a raw forecast file
#' @param verbose boolean Output warnings of specific violations
#' @return boolean
#'
forecast_is_valid <- function(forecast_df, verbose = FALSE){
  is_valid <- TRUE
  violations <- c()

  valid_columns <- c("date", "forecastmonth", "forecastyear", "newmoonnumber", 
                     "model", "currency", "level", "species", "estimate",
                     "LowerPI", "UpperPI", "fit_start_newmoon",
                     "fit_end_newmoon", "initial_newmoon")
  valid_currencies <- c("abundance", "richness", "biomass", "energy")
  valid_levels <- paste("Plot", 1:24, " ", sep = "")
  valid_levels <- c("All", "Controls", "FullExclosure", "KratExclosure", 
                    valid_levels)
  valid_species <- c("total", "BA", "DM", "DO", "DS", "OL", "OT", "PB", "PE",
                     "PF", "PH", "PI", "PL", "PM", "PP", "RF", "RM", "RO", 
                     "SF", "SH", "SO", "NA")

  if(!(all(colnames(forecast_df) %in% valid_columns) & 
       all(valid_columns %in% colnames(forecast_df)))){
    if(verbose){
      print("Forecast file column names invalid")
    }
    return(FALSE)
  }

  # TODO: Account for dates that are formatted correctly but potentially
  #  many years off.
  forecast_df$date <- as.Date(forecast_df$date, "%Y-%m-%d")
  if(any(is.na(forecast_df$date))){
    is_valid <- FALSE
    violations <- c("date", violations)
  }
  if(!all(unique(forecast_df$currency) %in% valid_currencies)){
    is_valid <- FALSE
    violations <- c("currency", violations)
  }
  if(!all(unique(forecast_df$level) %in% valid_levels)){
    is_valid <- FALSE
    violations <- c("level", violations)
  }
  if(!all(unique(forecast_df$species) %in% valid_species)){
    is_valid <- FALSE
    violations <- c("species", violations) 
  }
  if(any(is.na(forecast_df$estimate))){ 
    is_valid <- FALSE
    violations <- c("NA esimates", violations) 
  }
  if(any(is.na(forecast_df$LowerPI))){ 
    is_valid <- FALSE
    violations <- c("NA LowerPI", violations)
  }
  if(any(is.na(forecast_df$UpperPI))){
    is_valid <- FALSE 
    violations <- c("NA UpperPI", violations)
  }
  if(!is.integer(forecast_df$fit_start_newmoon)){ 
    is_valid <- FALSE
    violations <- c("fit_start_newmoon not int")
  }
  if(!is.integer(forecast_df$fit_end_newmoon)){ 
    is_valid <- FALSE
    violations <- c("fit_end_newmoon not int")
  }
  if(!is.integer(forecast_df$initial_newmoon)){ 
    is_valid <- FALSE
    violations <- c("initial_newmoon not int")
  }
  if(any(is.na(forecast_df$fit_start_newmoon))){ 
    is_valid <- FALSE
    violations <- c("fit_start_newmoon contains NA")
  }
  if(any(is.na(forecast_df$fit_end_newmoon))){ 
    is_valid <- FALSE
    violations <- c("fit_end_newmoon contains NA")
  }
  if(any(is.na(forecast_df$initial_newmoon))){ 
    is_valid <- FALSE
    violations <- c("initial_newmoon contains NA")
  }
  
  if(verbose & length(violations) > 0){
    print(paste("Forecast validation failed: ", violations), sep = "")
  }
  return(is_valid)
}


#' Collect all separate forecasts file into a single dataframe.
#'
#' The base folder can include subfolders.
#' Will only include files which pass validation. Will issue
#' warnings if a file isn't valid. If verbose is True it will
#' print specifics about the validation.
#'
#' @param forecast_folder str Base folder holding all forecast files
#' @param verbose bool Output info on file violations
#' @return dataframe combined forecasts
#'
compile_forecasts <- function(forecast_folder = './predictions', 
                              verbose = FALSE, use_hindcasts = FALSE){
  if (use_hindcasts){
    search_string <- "hindcast"
  } else{
    search_string <- "forecast"
  }
  
  fcast_filenames <- list.files(forecast_folder, pattern = search_string, 
                          full.names = TRUE, recursive = TRUE)
  fcast_filenames <- fcast_filenames[!grepl("model_aic", fcast_filenames)]
  all_forecasts <- data.frame()

  for(this_forecast_file in fcast_filenames){
    this_forecast_data <- try(read.csv(this_forecast_file, na.strings = "", 
                                stringsAsFactors = FALSE))
    if(class(this_forecast_data) %in% "try-error"){
      if(verbose){
        print(paste("File not readable: ", this_forecast_file, sep = ""))
      }else {
        warning(paste("File not readable: ", this_forecast_file, sep = ""))
      }
      next
    }

    if (verbose){
      print(paste("Testing file ", this_forecast_file, sep = ""))
    }
    if (forecast_is_valid(this_forecast_data, verbose = verbose)){
      if(verbose) {
        print(paste("File format is valid: ", this_forecast_file, sep = ""))
        print("-------")
      }
      all_forecasts <- all_forecasts %>%
                         dplyr::bind_rows(this_forecast_data)
    } else {
      if (verbose){
        print(paste("File format not valid: ", this_forecast_file, sep = ""))
        print("-------")
      } else{
        warning(paste("File format not valid: ", this_forecast_file, sep=""))
      }
    }
  }
  all_forecasts$date <- as.Date(all_forecasts$date)
  return(all_forecasts)
}

#' Visualize a time-series forecast
#'
#' @details Plots the observed time-series and the 1-step forecasts within it
#'   Plots the forecast time-series along with the prediction interval for 
#'   future observations
#' @param obs_data is a data.frame (observed data)
#' @param obs_date_col_name is a string: name of the date column from obs_data
#' @param obs_val_col_name is a string: name of the column of the value being
#'   forecast
#' @param for_data is a data.frame (forecast data)
#' @param for_date_col_name is a string: name of the date column from for_data
#' @param for_val_col_name is a string: name of the column of value being 
#'   forecast, from for_data
#' @param for_model_name is a string: name of the model to be used from model 
#'   column in for_data
#' @param for_lowerpi_col_name is a string: name of the column of the lower 
#'   confidence interval from for_data
#' @param for_upperpi_col_name is a string: name of the column of the upper 
#'   confidence interval from for_data
#' @param start_newmoon is numeric: first new moon number to be plotted
#' @param ylabel is a string: title for y-axis
#'
forecast_viz <- function(obs_data, obs_date_col_name, obs_val_col_name, 
                         for_data, for_date_col_name, for_val_col_name,
                         for_model_name, for_lowerpi_col_name, 
                         for_upperpi_col_name, start_newmoon, ylabel){

  for_data_sub <- dplyr::filter(for_data, species == obs_val_col_name, 
                    model == for_model_name)
  obs_data_sub <- dplyr::filter(obs_data, newmoonnumber >= start_newmoon)

  ggplot2::ggplot(obs_data_sub, ggplot2::aes_string(x = obs_date_col_name)) +
  ggplot2::geom_ribbon(data = for_data_sub, 
    mapping = ggplot2::aes_string(x = for_date_col_name, 
                ymin = for_lowerpi_col_name, ymax = for_upperpi_col_name), 
    fill = "lightblue") +
  ggplot2::geom_line(ggplot2::aes_string(y = obs_val_col_name)) +
  ggplot2::geom_line(data = for_data_sub, 
    mapping = ggplot2::aes_string(x = for_date_col_name, 
                y = for_val_col_name), 
    color = "blue") +
  ggplot2::labs(x = "", y = ylabel)
}
