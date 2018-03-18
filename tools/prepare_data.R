"%>%" <- magrittr::"%>%"
source("tools/data_tools.R")

filename_suffix <- "forecasts"

# The date this forecast is run. Always today's date.
forecast_date <- Sys.Date()

# Gather the data
portalr::download_observations(release_only = FALSE)

moons <- get_moon_data()
curr_moons <- moons %>%
              dplyr::select(newmoonnumber, newmoondate, period, censusdate)
curr_moons$newmoondate <- as.Date(as.character(curr_moons$newmoondate))
future_moons <- portalr::get_future_moons(moons)
total_moons <- rbind(curr_moons, future_moons)

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

# Set up the beginning and end of the forecast timeperiods
#  this acknowledges that the rodents might not have been sampled on the most
#  recent newmoons, even though we now have covariate values for them  
#
# prev = previous (i.e. the most recent)

which_prev_newmoon <- max(which(total_moons$newmoondate < forecast_date))
prev_newmoon <- moons$newmoonnumber[which_prev_newmoon]

prev_rodent_period_all <- tail(rodent_data$all, 1)$period
prev_rodent_period_control <- tail(rodent_data$control, 1)$period
prev_rodent_period <- max(prev_rodent_period_all, prev_rodent_period_control)
which_prev_rodent_period <- which(moons$period == prev_rodent_period)
prev_rodent_newmoon <- moons$newmoonnumber[which_prev_rodent_period]

prev_covar_nm_weather <- tail(weather_data, 1)$newmoonnumber
prev_covar_nm_ndvi <- tail(ndvi_data, 1)$newmoonnumber
prev_covar_newmoon <- min(prev_covar_nm_weather, prev_covar_nm_ndvi)

first_fcast_covar_newmoon <- prev_covar_newmoon + 1
first_fcast_rodent_newmoon <- prev_rodent_newmoon + 1
last_fcast_newmoon <- prev_newmoon + 12

rodent_fcast_newmoons <- first_fcast_rodent_newmoon:last_fcast_newmoon
which_r_nms <- which(total_moons$newmoonnumber %in% rodent_fcast_newmoons)
rodent_nm_dates <- total_moons$newmoondate[which_r_nms]
rodent_fcast_months <- as.numeric(format(rodent_nm_dates, "%m"))
rodent_fcast_years <- as.numeric(format(rodent_nm_dates, "%Y"))

covar_fcast_newmoons <- first_fcast_covar_newmoon:last_fcast_newmoon
which_c_nms <- which(total_moons$newmoonnumber %in% covar_fcast_newmoons)
covar_nm_dates <- total_moons$newmoondate[which_c_nms]
covar_fcast_months <- as.numeric(format(covar_nm_dates, "%m"))
covar_fcast_years <- as.numeric(format(covar_nm_dates, "%Y"))

# Write data files
write.csv(rodent_data$all, "data/rodent_all.csv", row.names = FALSE)
write.csv(rodent_data$controls, "data/rodent_controls.csv", row.names = FALSE)
write.csv(weather_data, "data/weather_data.csv", row.names = FALSE)
write.csv(moons, "data/moon_data.csv", row.names = FALSE)
write.csv(ndvi_data, "data/ndvi_data.csv", row.names = FALSE)

# Write YAML
yaml_list <-  list(filename_suffix = filename_suffix, 
                forecast_date = as.character(forecast_date), 
                covariate_forecast_newmoons = covar_fcast_newmoons, 
                covariate_forecast_months = covar_fcast_months, 
                covariate_forecast_years = covar_fcast_years,
                rodent_forecast_newmoons = rodent_fcast_newmoons, 
                rodent_forecast_months = rodent_fcast_months, 
                rodent_forecast_years = rodent_fcast_years)
writeLines(yaml::as.yaml(yaml_list), con = "data/model_metadata.yaml")
