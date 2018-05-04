"%>%" <- magrittr::"%>%"
"!!!" <- rlang::"!!!"
"!!" <- rlang::"!!"

source("tools/data_tools.R")

filename_suffix <- "forecasts"
forecast_date <- Sys.Date()
portalr::download_observations()

moons <- get_moon_data()
total_moons <- add_future_moons(moons, forecast_date)
rodent_data <- get_rodent_data(moons, forecast_date)
weather_data <- get_weather_data()
ndvi_data <- portalr::ndvi("newmoon", fill = TRUE)

write.csv(rodent_data$all, "data/rodent_all.csv", row.names = FALSE)
write.csv(rodent_data$controls, "data/rodent_controls.csv" ,row.names = FALSE)
write.csv(weather_data, "data/weather_data.csv", row.names = FALSE)
write.csv(moons, "data/moon_data.csv", row.names = FALSE)
write.csv(ndvi_data, "data/ndvi_data.csv", row.names = FALSE)

model_metadata <- prep_metadata(rodent_data, total_moons)
writeLines(yaml::as.yaml(model_metadata), con = "data/model_metadata.yaml")