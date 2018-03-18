#
# I am actively working to update this script based on the new covariate 
#  forecasting. as of now, the code is not functional with the updated model
#  scripts and metadata structure
#

source("tools/forecast_tools.R")
source("tools/data_tools.R")
"%>%" <- magrittr::"%>%"

filename_suffix <- "hindcasts"

# The date this hindcast is run. Always today's date.
forecast_date <- Sys.Date()

# Hindcast will set the time based on these newmoon numbers. For each one,
# a hindcast will be made which pretends that sampling period had just 
# happened. #403 to 490 is Jan,2010 - Jan,2017. 
initial_time_newmoons <- 490:403

# Get the most recent data loaded into the data folder
portalr::download_observations(release_only = FALSE)
moons <- get_moon_data()
rodent_data <- get_rodent_data(moons, forecast_date)
cols <- c("mintemp", "maxtemp", "meantemp", "precipitation", "newmoonnumber")
weather_data <- portalr::weather("newmoon", fill = TRUE) %>% 
                dplyr::ungroup() %>%
                dplyr::select(cols)
incompletes <- which(is.na(weather_data$newmoonnumber))
if (length(incompletes) > 0){
  weather_data <- weather_data[-incompletes, ]
}
ndvi_data <- portalr::ndvi("newmoon", fill = TRUE)

# Write data files
write.csv(rodent_data$all, "data/rodent_all.csv", row.names = FALSE)
write.csv(rodent_data$controls, "data/rodent_controls.csv", row.names = FALSE)
write.csv(weather_data, "data/weather_data.csv", row.names = FALSE)
write.csv(moons, "data/moon_data.csv", row.names = FALSE)
write.csv(ndvi_data, "data/ndvi_data.csv", row.names = FALSE)

trap_refpath <- "PortalData/Rodents/Portal_rodent_trapping.csv"
trap_fpath <- portalr::FullPath(trap_refpath, "~")
trappings <- read.csv(trap_fpath)
incomplete_samples <- portalr::find_incomplete_censuses(trappings)

failed_newmoons <- c()

for(this_newmoon in initial_time_newmoons){

  moons <- get_moon_data()
  
  this_censusdate <- moons %>%
                     dplyr::filter(newmoonnumber == this_newmoon) %>%
                     dplyr::pull(censusdate) %>%
                     as.character() %>%
                     as.Date()
  
  this_period <- moons %>%
                 dplyr::filter(newmoonnumber == this_newmoon) %>%
                 dplyr::pull(period)
  
  # Don't do hindcasting from newmoons that had incomplete samplings
  # or were not sampled at all
  this_sample_incomplete <- this_period %in% incomplete_samples$period 
  if (this_sample_incomplete| is.na(this_censusdate)){
    next
  }
  
  # Set up the data for hindcasting
  moons <- get_moon_data() %>%
           dplyr::filter(newmoonnumber<=this_newmoon)
  curr_moons <- moons %>%
                dplyr::select(newmoonnumber, newmoondate, period, censusdate)
  curr_moons$newmoondate <- as.Date(as.character(curr_moons$newmoondate))
  future_moons <- portalr::get_future_moons(moons)
  total_moons <- rbind(curr_moons, future_moons)
  
  # Beginning and end of the forecast timeperiod
  most_recent_newmoon <- moons$newmoonnumber[which.max(moons$period)]
  most_recent_newmoon_date <- as.Date(moons$newmoondate[which.max(moons$period)])
  first_fcast_newmoon <- most_recent_newmoon + 1
  last_fcast_newmoon <- first_fcast_newmoon + 11
  fcast_nms <- first_fcast_newmoon:last_fcast_newmoon
  which_in <- which(future_moons$newmoonnumber %in% fcast_nms)
  fcast_dates <- future_moons$newmoondate[which_in]
  fcast_months <- format(fcast_dates, "%m") 
  fcast_years <- format(fcast_dates, "%Y") 
  
  all <- read.csv("data/rodent_all.csv") %>%
         dplyr::filter(newmoonnumber <= this_newmoon)
  
  controls <- read.csv("data/rodent_controls.csv") %>%
              dplyr::filter(newmoonnumber <= this_newmoon)
  
  weather_data <- read.csv("data/weather_data.csv") %>%
                 dplyr::filter(newmoonnumber <= this_newmoon)
  
  # Update files in tools directory to use in this specific hindcast
  # Write data files
  write.csv(all, "data/rodent_all.csv", row.names = FALSE)
  write.csv(controls, "data/rodent_controls.csv", row.names = FALSE)
  write.csv(weather_data, "data/weather_data.csv", row.names = FALSE)
  write.csv(ndvi_data, "data/ndvi_data.csv", row.names = FALSE)
  
  # Write YAML
  writeLines(
    yaml::as.yaml(list(filename_suffix = filename_suffix,
      forecast_date = as.character(forecast_date), 
      forecast_newmoons = forecast_newmoons, 
      forecast_months = forecast_months, 
      forecast_years = forecast_years)), con = "data/model_metadata.yaml")
  
  # Run all models
  cat("Running models", "\n")
  dir.create("tmp")

  print("######################################")
  print("######################################")
  print(paste0("Running models for initial newmoon: ", this_newmoon))
  print("######################################")
  print("######################################")

  model_outcome <- try(sapply(list.files("models", full.names = TRUE), source))
  if (class(model_outcome) == "try-error"){
    print("######################################")
    print(paste0("Failed doing hindcasts for newmoon: ", this_newmoon)) 
    print("######################################")
    failed_newmoons <- c(failed_newmoons, this_newmoon)
  }
 
  # Compile all hindcasts into one file
  allhindcasts <- forecastall(forecast_date, filename_suffix)
  unlink("tmp/*")
}

print_msg <- paste0("Hindcasting complete. Failed on ", 
               length(failed_newmoons)," initial newmoons: ")
print(print_msg)
print(failed_newmoons)
