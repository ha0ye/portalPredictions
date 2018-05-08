#' Combine all new forecasts (from the tmp directory), add ensembles
#' 
#' @param forecast_date
#' @param filename_suffix
#' @return list(forecasts,all_model_aic)
#' @example forecastall('forecasts')
#'
forecastall <- function(forecast_date, filename_suffix = "forecasts"){
  
  file_ptn <- paste(filename_suffix, ".csv", sep = "")
  files <- list.files("tmp", pattern = file_ptn, full.names = TRUE)
  col_class <- c("Date", "integer", "integer", "integer", "character", 
                 "character", "character", "character", "numeric", "numeric",
                 "numeric", "integer", "integer", "integer")
  fcasts <- do.call(rbind, 
            lapply(files, read.csv, na.strings = "", colClasses  = col_class))

  file_ptn <- paste(filename_suffix, "_model_aic.csv", sep = "")
  files <- list.files("tmp", pattern = file_ptn, full.names = TRUE)

  aics <- do.call(rbind,
            lapply(files, read.csv, na.strings = ""))
  
  fcast_date <- as.character(forecast_date)
  fcast_fname <- paste(fcast_date, filename_suffix, ".csv", sep = "")
  forecast_filename <- file.path("predictions", fcast_fname)
  aic_fname <- paste(fcast_date, filename_suffix, "_model_aic.csv", sep = "")
  model_aic_filename <- file.path("predictions", aic_fname)
  append_csv(fcasts, forecast_filename)
  append_csv(aics, model_aic_filename)
  
  ensemble <- make_ensemble(fcasts) %>% 
                subset(select = colnames(fcasts))
  append_csv(ensemble, forecast_filename)
  
  return(list(forecasts = fcasts, all_model_aic = aics))
}

#' Appending a csv without re-writing the header.
#' @param df data table to be written out
#' @param filename filename of existing csv to be appended
#' @return 
#' @export
#'
append_csv <- function(df, filename){
  write.table(df, filename, sep = ",", row.names = FALSE, 
    col.names = !file.exists(filename), append = file.exists(filename))
}

#' calculate akaike weights
#' @param forecast_folder folder where the forecast files are
#' @return model weights
#' @export
#'
compile_aic_weights <- function(forecast_folder = "./predictions"){
  aic_files <- list.files(forecast_folder, full.names = TRUE, 
                          recursive = TRUE)
  aic_files <- aic_files[grepl("model_aic", aic_files)]

  aics <- purrr::map(aic_files, 
            ~read.csv(.x, na.strings = "", stringsAsFactors = FALSE)) %>% 
          dplyr::bind_rows()
 
  grps <- rlang::quos(date, currency, level, species, fit_start_newmoon, 
            fit_end_newmoon, initial_newmoon)
 
  wts <- aics %>%
    dplyr::group_by(!!!grps) %>%
    dplyr::mutate(delta_aic = aic - min(aic), 
              weight = exp(-0.5 * delta_aic) / sum(exp(-0.5*delta_aic))) %>%
    dplyr::ungroup()

  return(wts)
}

#' Create the ensemble model from all other forecasts
#' 
#' @description Uses the weighted mean and weighted sample variance
#'   Mean is the weighted mean of all model means. Variance is the weighted 
#'   mean of all model variances + the variances of the weighted mean using
#'   the unbiased estimate of sample variance. See
#'   https://github.com/weecology/portalPredictions/pull/65
#'   We only store the prediction interval for models, so backcalculate 
#'   individual model variance assuming the same CI_level throughout. Assert
#'   that the summed weight of all the model ensembles is 1, as that's what
#'   the above variance estimates assume. Rounded to account for precision
#'   errors. Summed weights can also be NA if there are not weights availble
#'   for that ensemble. 
#' @param all_forecasts alll forecasts
#' @param models_to_use models to use
#' @param CI_level confidence interval level
#' @return ensemble
#' @export
#'
make_ensemble <- function(all_forecasts, models_to_use = NA, CI_level = 0.9){

  weights <- compile_aic_weights()
  weights$date <- as.Date(weights$date)
  CI_quantile <- qnorm((1 - CI_level) / 2, lower.tail = FALSE)

  leftj <- c("date", "model", "currency", "level", "species", 
             "fit_start_newmoon", "fit_end_newmoon", "initial_newmoon")
  grp <- rlang::quos(date, newmoonnumber, forecastmonth, forecastyear,level, 
           currency, species, fit_start_newmoon, fit_end_newmoon, 
           initial_newmoon)
  mod_var <- rlang::quo(((UpperPI - estimate)/CI_quantile) ^ 2)
  est <- rlang::quo(sum(estimate * weight))
  wtss <- rlang::quo(sum(weight * (estimate - ensemble_estimate)^2))
  ens_var <- rlang::quo(sum(model_var * weight) + 
                        weighted_ss / (n() * sum(weight) - 1))
  weighted_estimates <- all_forecasts %>%
                        dplyr::mutate(model_var = !!mod_var) %>%
                        dplyr::left_join(weights, by = leftj) %>%
                        dplyr::group_by(!!!grp) %>%
                        dplyr::summarise(ensemble_estimate = !!est, 
                                         weighted_ss = !!wtss ,
                                         ensemble_var = !!ens_var,
                                         sum_weight = sum(weight)) %>% 
                        dplyr::ungroup() 
     
  check_sum <- round(weighted_estimates$sum_weight, 10) == 1 
  check_na <- is.na(weighted_estimates$sum_weight)       
  if(!all(check_sum | check_na)){ 
    stop("Summed weights do not equal 1")
  }

  PI_l <- rlang::quo(ensemble_estimate - (sqrt(ensemble_var) * CI_quantile))
  PI_U <- rlang::quo(ensemble_estimate + (sqrt(ensemble_var) * CI_quantile))
  ensemble <- weighted_estimates %>%
              dplyr::mutate(LowerPI = !!PI_l, UpperPI = !!PI_U) %>%
              dplyr::mutate(LowerPI = ifelse(LowerPI<0, 0, LowerPI)) %>%
              dplyr::rename(estimate = ensemble_estimate) %>%
              dplyr::select(-ensemble_var, -weighted_ss, -sum_weight)
  
  ensemble$model <- "Ensemble"
  return(ensemble)
}


#' Get sp predictions
#' @param data data
#' @param lvl lvl
#' @param lead_time lead_time
#'
#' @export
#'
get_sp_predicts <- function(data, lvl, lead_time) {
  data <- data %>% 
          transform(
            forecast_date = 
            zoo::as.yearmon(paste(forecastmonth, "/", forecastyear, sep = ""),
            format = "%m/%Y")) %>% 
          transform(date = as.Date(date, "%Y-%m-%d"))


  data1 <- dplyr::filter(data, level == lvl,
                 date == max(as.Date(date)))
  target_moon <-min(data1$newmoonnumber) + (lead_time - 1)
  data2 <- dplyr::filter(data1, newmoonnumber == target_moon)
}

#' Plot species forecast
#'
#' @description produces the second plot on the Species-level Forecast on the
#'   Current Forecast page on the website
#'  
#' @param data
#' @param title main title for plot
#' @return sp_predict is a plot object -- plot(sp_predict) displays it
#' 
plot_species_forecast <- function(data, title){

  url_base <- "https://raw.githubusercontent.com/weecology/PortalData/"
  url_suff <- "master/Rodents/moon_dates.csv"
  url <- paste(url_base, url_suff, sep = "")
  newmoons_table <- read.csv(text = RCurl::getURL(url))
  target_moon <- unique(data$newmoonnumber)
  period_code <- newmoons_table %>%
                 dplyr::filter(newmoonnumber == target_moon) %>%
                 dplyr::select(period) %>%
                 as.integer()

  url_suff <- "master/Rodents/Portal_rodent_species.csv"
  url <- paste(url_base, url_suff, sep = "")
  species_table <- read.csv(text = RCurl::getURL(url),  
                            stringsAsFactors = FALSE, na.strings = "")
  species_names <- species_table %>% 
                   dplyr::select("speciescode", "scientificname") %>% 
                   rbind(c("total", "total")) %>%
                   merge(data[ , c("species", "estimate")], 
                         by.x = "speciescode", by.y= "species")
  
  sp_predict <- ggplot2::ggplot(data,
                                ggplot2::aes(
                                  x = estimate,
                                  y = reorder(species, estimate),
                                  xmin = LowerPI,
                                  xmax = UpperPI
                                )) +
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

#' Compare forecasts to observations over different lead times.
#' 
#' @description Error can be any function. The level, species, and currency 
#'   columns from observations and forecasts must have matching values. Note
#'   this gives an average error value over many forecast iterations. Will
#'   only return values where there are matching comparison columns 
#'   (currency, level, species)
#'
#' @param observations dataframe Has the columns newmoonnumber, currency, 
#'   level, species, actual
#' @param forecasts dataframe passes the forecast validity check. Must have 
#'   matching values in the comparison columns
#' @param error_metric chr either "RMSE" for root mean squared error or 
#'   "coverage" for the coverage of the prediction intervals
#' @param ci_value int The value of the forecast confidence interval to scale
#'   PI values for the likelihood metric
#' @return data.frame Data.frame with the columns model, error, lead_time, 
#'   level, species, currency
#'


# working on this function right now. need to make sure it works
# as i go along tho, so next step is to do that then jump back into it


calculate_forecast_error <- function(observations, forecasts, 
                                     error_metric = "RMSE", CI_level = 0.9){

  #The tibble datatype output from dplyr causes issues here
  observations <- as.data.frame(observations)
  forecasts <- as.data.frame(forecasts)

  if(!forecast_is_valid(forecasts)){
    stop("Forecast dataframe not valid")
  }

  needed_cols <- c("newmoonnumber", "currency", "level", "species", "actual")
  if(!all(needed_cols %in% colnames(observations))){
    stop("observation data.frame does not have valid column names")
  }

  #At least 1 matching value must be in each of these columns in the 
  # observations and forecasts
  #TODO: Ensure matching rows in all 3 columns at once instead of just 
  # one at a time.

  column_check <- c()
  for(coli in c("currency", "level", "species")){
    if(!any(unique(observations[ , coli]) %in% unique(forecasts[ , coli]))){
      column_check <- c(column_check, column)
    }
  }
  if(length(column_check) >0){
    msg <- paste("Columns do not match: ", column_check, collapse = " ")
    stop(msg)
  }

  #Summarize to mean error by lead time. Lead time is number of new moons
  # ahead of when the forecast was made. This assumes a forecast was made with
  # only the data available prior to the first NewMoonDate in the series.
  #TODO: Make the lead time the actual days or weeks once more frequent 
  # forecasts are being made( see #37)

  joins <- c("newmoonnumber", "currency", "level", "species")

  forecasts <- forecasts %>%
               dplyr::mutate(lead_time = newmoonnumber - initial_newmoon) %>%
               dplyr::inner_join(observations, by = joins)

  if(error_metric == "RMSE"){
    grouping <- rlang::quos(model, currency, level, species, lead_time)
    comparisons <- forecasts %>%
                   dplyr::mutate(error_value = (estimate - actual)^2) %>%
                   dplyr::group_by(!!!grouping) %>%
                   dplyr::summarize(error_value = sqrt(mean(error_value))) %>%
                   dplyr::ungroup() %>%
                   dplyr::mutate(error_metric = "RMSE")

  } else if(error_metric == "coverage") {
    grouping <- rlang::quos(model, currency, level, species, lead_time, 
                            error_metric)
    comparisons <- forecasts %>%
                   dplyr::mutate(
                     within_prediction_interval = actual >= LowerPI & 
                       actual <= UpperPI, 
                     error_metric = "coverage") %>%
                   dplyr::group_by(!!!grouping) %>%
                   dplyr::summarize(
                     error_value = mean(within_prediction_interval)) %>%
                   dplyr::ungroup() %>%
                   dplyr::mutate(error_metric = "coverage")
    
  } else if(error_metric == "deviance") {
    stop("Deviance not implimented  yet")
  } else {
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
#' @return plot
#' @export
#'
plot_lead_time_errors <- function(error_df, level, species, currency, 
                                  error_metric){
  plot_title <- paste0("Level: ", level, ", Species: ", species, 
                       ", Currency: ", currency)

  graph <- ggplot2::ggplot(error_df, 
                           ggplot2::aes(x = lead_time, y = error, 
                                        group = model, color = model)) +
           ggplot2::geom_point()+
           ggplot2::geom_line() +
           ggplot2::labs(y = error_metric, x = "Lead Time (New Moons)", 
                title = plot_title)
  plot(graph)
}

#' Ensure that a forecast file is in the correct format
#'
#' @description Tools for working with forecast data expect a certain format.
#' This ensures a forecast file meets those formats. All column
#' and variable names are case sensitive. Column names should match exactly, 
#' case and everything, no more, no less. For specification see:
#' https://github.com/weecology/portalPredictions/wiki/forecast-file-format
#' Dates must be in the format YYYY-MM-DD
#'
#' @param forecast_df dataframe A dataframe read from a raw forecast file
#' @param verbose boolean Output warnings of specific violations
#' @return boolean
#' @export
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

  all_ok <- all(colnames(forecast_df) %in% valid_columns)
  all_incl <-  all(valid_columns %in% colnames(forecast_df))
  if(!(all_ok & all_incl)){
    if(verbose){
      print("Forecast file column names invalid")
    }
    return(FALSE)
  }

  #TODO: Account for dates that are formatted correctly but potentially many 
  # years off.
  forecast_df$date <- base::as.Date(forecast_df$date, "%Y-%m-%d")
  if(any(is.na(forecast_df$date))){
    is_valid <- FALSE
    violations <- c("date", violations)
  }

  #All of these must be ones listed in the forecast format wiki.
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

  #Estimates and PI's cannot have NA values
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

  #All the newmoon columns should be whole numbers with nothing missing
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
    print(paste("Forecast validation failed: ", violations), sep="")
  }
  return(is_valid)
}


#' Collect all separate forecasts file into a single dataframe.
#'
#' @description The base folder can include subfolders. Will only include
#'   files which pass validation. Will issue warnings if a file isn't valid.
#'   If verbose is True it will print specifics about the validation.
#'
#' @param forecast_folder str Base folder holding all forecast files
#' @param verbose bool Output info on file violations
#' @return dataframe combined forecasts
#'
#' @export
#'
compile_forecasts <- function(forecast_folder = "./predictions", 
                              verbose = FALSE, use_hindcasts = FALSE){
  if(use_hindcasts){
    search_string <- "hindcast"
  } else {
    search_string <- "forecast"
  }
  
  forecast_filenames <- list.files(forecast_folder, pattern = search_string, 
                                   full.names = TRUE, recursive = TRUE)
  forecast_filenames <- forecast_filenames[
                          !grepl("model_aic", forecast_filenames)]
  all_forecasts <- data.frame()

  for(this_forecast_file in forecast_filenames){
    this_forecast_data <- try(
                            read.csv(this_forecast_file, na.strings = "", 
                            stringsAsFactors = FALSE))
    if(class(this_forecast_data) %in% "try-error"){
      if(verbose){
        print(paste("File not readable: ", this_forecast_file, sep=""))
      } else {
        warning(paste("File not readable: ", this_forecast_file, sep=""))
      }
      next
    }

    if(verbose){
      print(paste("Testing file ", this_forecast_file, sep=""))
    }
    if(forecast_is_valid(this_forecast_data, verbose=verbose)){
      if(verbose) {
        print(paste("File format is valid: ", this_forecast_file, sep=""))
        print("-------")
      }
      all_forecasts <- dplyr::bind_rows(all_forecasts, this_forecast_data)
    } else {
      if(verbose){
        print(paste("File format not valid: ", this_forecast_file, sep=""))
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
#' @description Plots the observed time-series and the 1-step forecasts within
#'   it. Plots the forecast time-series along with the prediction interval for
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
#' @export
#'
forecast_viz <- function(obs_data, obs_date_col_name, obs_val_col_name, 
                         for_data, for_date_col_name, for_val_col_name, 
                         for_model_name, for_lowerpi_col_name, 
                         for_upperpi_col_name, start_newmoon, ylabel){

  for_data_sub <- dplyr::filter(for_data, species == obs_val_col_name, 
                                model == for_model_name)
  obs_data_sub <- dplyr::filter(obs_data, newmoonnumber >= start_newmoon)
  max_obs <- max(obs_data_sub$newmoonnumber)
  min_for <- min(for_data_sub$newmoonnumber)
  which_max_obs <- which(obs_data_sub$newmoonnumber == max_obs)
  which_min_for <- which(for_data_sub$newmoonnumber == min_for)
  last_obs <- as.numeric(obs_data_sub[which_max_obs, obs_val_col_name])
  first_for <- for_data_sub[which_min_for, for_val_col_name]
  last_obs_date <- obs_data_sub[which_max_obs, obs_date_col_name]
  first_for_date <- for_data_sub[which_min_for, for_date_col_name]
  conx <- c(data.frame(last_obs_date)[1, 1], first_for_date)
  cony <- c(last_obs, first_for)
  conxy <- data.frame(conx, cony)

  ggplot2::ggplot(obs_data_sub, ggplot2::aes_string(x = obs_date_col_name)) +
                  ggplot2::geom_ribbon(data = for_data_sub, 
                  mapping = ggplot2::aes_string(x = for_date_col_name, 
                                                ymin = for_lowerpi_col_name,
                                                ymax = for_upperpi_col_name), 
                                                fill = "lightblue") +
  ggplot2::geom_line(ggplot2::aes_string(y = obs_val_col_name)) +
  ggplot2::geom_line(data = for_data_sub, 
                     mapping = ggplot2::aes_string(x = for_date_col_name, 
                                                   y = for_val_col_name), 
                     color = "blue") +
  ggplot2::labs(x = "", y = ylabel) +
  ggplot2::geom_line(data = conxy, lty = 3,
                     mapping = ggplot2::aes_string(x = conx, y = cony))
}
