#' Get data for newmoonnumbers and related trapping period codes
#' @return newmoons table
#' @examples
#' get_moon_data()
#' 
get_moon_data <- function(){
  file_path <- portalr::FullPath('PortalData/Rodents/moon_dates.csv', '~')
  moons <- read.csv(file_path, header = T)
  moons$year <- lubridate::year(moons$newmoondate)
  moons$month <- lubridate::month(moons$newmoondate)
  return(moons)
}

#' Get rodent data, tailored for forecasting (all plots and controls only)
#' @details Training data currently being in Jan 2995 (period 203/newmoon 217)
#' @param moons current newmoonnumber table
#' @param forecast_date date the forecast is run
#' 
#' @return a list of two dataframes, all plots and control plots
#'
get_rodent_data <- function(moons, forecast_date){
  historic_start_period <- 203
  historic_start_newmoon <- 217
  
  # Control plots
  controls <- portalr::abundance(level = "Treatment", type = "Rodents",
                length = "Longterm", incomplete = FALSE)
  # Drop PI
  controls <- controls[ , -which(colnames(controls) == "PI")]

  # The total rodent count in each treatment
  controls$total <- rowSums(controls[ , -(1:2)])

  # Drop non-control treatments and add in newmoonnumber
  controls <- controls %>%
                dplyr::filter(treatment == 'control') %>%
                dplyr::select(-treatment) %>%
                dplyr::inner_join(moons, by = c("period" = "period")) %>%
                subset(newmoonnumber >= historic_start_newmoon) %>%
                dplyr::select(-newmoondate, -censusdate)
  
  # All plots
  all <- portalr::abundance(level = "Site", type = "Rodents", length = "all", 
           incomplete = FALSE)
  # Drop PI
  all <- all[ , -which(colnames(all) == "PI")]

  # The total rodent count across the entire site
  all$total <- rowSums(all[ , -(1)])
  all <- all %>% 
           dplyr::inner_join(moons,by = c("period"="period")) %>%
           subset(period >= historic_start_period) %>%
           dplyr::select(-newmoondate, -censusdate)
    
  rodent_data <- list()
  rodent_data$controls <- controls
  rodent_data$all <- all
  return(rodent_data)
}
