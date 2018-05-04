#' Return normalized path for all operating systems
#' @param ReferencePath a path to join with current working directory
#' @param BasePath Current working directory else path given
#'
#' @return Path
#' @export
#' @examples
#' FullPath("PortalData/Rodents/Portal_rodent.csv")
#' FullPath("PortalData/Rodents/Portal_rodent.csv", "~")
#'
FullPath <- function(ReferencePath, BasePath = getwd()){
  BasePath <- normalizePath(BasePath)
  Path <- normalizePath(file.path(BasePath, ReferencePath), mustWork = FALSE)
  return(Path)
}

#' Get dates and data for newmoonnumbers and related trapping periods
#' @return newmoons table
#' @examples
#' get_moon_data()
#' 
get_moon_data <- function(){
  fpath <- FullPath("PortalData/Rodents/moon_dates.csv", "~")
  moons <- read.csv(fpath, header = TRUE)
  
  moons$year <- lubridate::year(moons$newmoondate)
  moons$month <- lubridate::month(moons$newmoondate)
  moons$newmoondate <- as.Date(as.character(moons$newmoondate))
  return(moons)
}

#' 
#' @param moons the moon data table
#' @param forecast_date the date of the forecast
#' @return moon data table with the appropriately extended dates 
#' @export
#'
add_future_moons <- function(moons, forecast_date){
  curr_moons <- moons %>%
                dplyr::select(newmoonnumber, newmoondate, period, censusdate)
  future_moons <- portalr::get_future_moons(moons)
  total_moons <- rbind(curr_moons, future_moons)

  not_future_moons <- which(future_moons$newmoondate < forecast_date)
  nnfm <- length(not_future_moons)
  if (nnfm > 0){
    add_nfm <- future_moons[not_future_moons, ]
    add_nfm$year <- as.numeric(format(add_nfm$newmoondate, "%Y"))
    add_nfm$month <- as.numeric(format(add_nfm$newmoondate, "%m"))
    moons <- rbind(moons, add_nfm)
    curr_moons <- moons %>% 
                  dplyr::select(newmoonnumber, newmoondate, period, censusdate)
    future_moons <- portalr::get_future_moons(moons)
    total_moons <- rbind(curr_moons, future_moons)
  }

  return(total_moons)
}

#' @return weather data table 
#' @export
#'
get_weather_data <- function(){
  cols <- c("mintemp", "maxtemp", "meantemp", "precipitation", "newmoonnumber")
  weather_data <- portalr::weather("newmoon", fill = TRUE) %>% 
                  dplyr::ungroup() %>%
                  dplyr::select(cols)
  incompletes <- which(is.na(weather_data$newmoonnumber))
  if (length(incompletes) > 0){
    weather_data <- weather_data[-incompletes, ]
  }
  return(weather_data)
}

#' Get rodent data, tailored for forecasting (all plots and controls only)
#' @param moons current newmoonnumber table
#' @param forecast_date date the forecast is run
#' 
#' @return a list of two dataframes, all plots and control plots
#' @examples
#' get_rodent_data(moons, forecast_date)
#'
get_rodent_data <- function(moons, forecast_date){

  # Corresponding to Jan 1995
  historic_start_period <- 203
  historic_start_newmoon <- 217
  
  # Control plot data
  controls <- portalr::abundance(clean = FALSE, level = "Treatment",
                                 type = "Rodents", length = "Longterm",
                                 min_plots = 24)
  # Drop PI
  controls <- controls[ , -which(colnames(controls) == "PI")]
  # The total rodent count in each treatment
  controls$total = rowSums(controls[,-(1:2)])
  # Drop non-control treatments and add in newmoonnumber
  controls <- controls %>%
              dplyr::filter(treatment == 'control') %>%
              dplyr::select(-treatment) %>%
              dplyr::inner_join(moons, by = c("period" = "period")) %>%
              subset(newmoonnumber >= historic_start_newmoon) %>%
              dplyr::select(-newmoondate, -censusdate)
  
  # All plot data
  all <- portalr::abundance(clean = FALSE, level = "Site", type = "Rodents",
                            length = "all", min_plots = 24)
  # Drop PI
  all <- all[ , -which(colnames(all) == "PI")]
  # The total rodent count across the entire site
  all$total <- rowSums(all[,-(1)])
  all <- all %>% 
         dplyr::inner_join(moons, by = c("period" = "period")) %>%
         subset(period >= historic_start_period) %>%
         dplyr::select(-newmoondate, -censusdate)
    
  rodent_data <- list()
  rodent_data$controls <- controls
  rodent_data$all <- all
  return(rodent_data)
}

#' Set up the model metadata, primarily the forecast timeperiods
#'
#' @param
#' @return model metadata as a list
#' @export
#'
prep_metadata <- function(rodents, moons){
 
  # prev = previous (i.e. the most recent)

  which_prev_newmoon <- max(which(moons$newmoondate < forecast_date))
  prev_newmoon <- moons$newmoonnumber[which_prev_newmoon]

  prev_rodent_pd_all <- tail(rodents$all, 1)$period
  prev_rodent_pd_control <- tail(rodents$control, 1)$period
  prev_rodent_pd <- max(prev_rodent_pd_all, prev_rodent_pd_control)
  which_prev_rodent_pd <- which(moons$period == prev_rodent_pd)
  prev_rodent_newmoon <- moons$newmoonnumber[which_prev_rodent_pd]

  prev_covar_nm_weather <- tail(weather_data, 1)$newmoonnumber
  prev_covar_nm_ndvi <- tail(ndvi_data, 1)$newmoonnumber
  prev_covar_newmoon <- min(prev_covar_nm_weather, prev_covar_nm_ndvi)

  first_fcast_covar_newmoon <- prev_covar_newmoon + 1
  first_fcast_rodent_newmoon <- prev_rodent_newmoon + 1
  last_fcast_newmoon <- prev_newmoon + 12

  rodent_fcast_newmoons <- first_fcast_rodent_newmoon:last_fcast_newmoon
  which_r_nms <- which(moons$newmoonnumber %in% rodent_fcast_newmoons)
  rodent_nm_dates <- moons$newmoondate[which_r_nms]
  rodent_fcast_months <- as.numeric(format(rodent_nm_dates, "%m"))
  rodent_fcast_years <- as.numeric(format(rodent_nm_dates, "%Y"))

  covar_fcast_newmoons <- first_fcast_covar_newmoon:last_fcast_newmoon
  which_c_nms <- which(moons$newmoonnumber %in% covar_fcast_newmoons)
  covar_nm_dates <- moons$newmoondate[which_c_nms]
  covar_fcast_months <- as.numeric(format(covar_nm_dates, "%m"))
  covar_fcast_years <- as.numeric(format(covar_nm_dates, "%Y"))

  out <-  list(filename_suffix = filename_suffix, 
               forecast_date = as.character(forecast_date), 
               covariate_forecast_newmoons = covar_fcast_newmoons, 
               covariate_forecast_months = covar_fcast_months, 
               covariate_forecast_years = covar_fcast_years,
               rodent_forecast_newmoons = rodent_fcast_newmoons, 
               rodent_forecast_months = rodent_fcast_months, 
               rodent_forecast_years = rodent_fcast_years)
  return(out)
}
