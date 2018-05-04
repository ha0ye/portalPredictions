source("tools/data_tools.R")
source("tools/forecast_tools.R")

model_metadata <- yaml::yaml.load_file("data/model_metadata.yaml")
forecast_date <- as.Date(model_metadata$forecast_date)

if (!(forecast_date == Sys.Date())){ 
  stop("Data not updated") 
}

cat("Running models", "\n")
dir.create("tmp")
sapply(list.files("models", full.names = TRUE), source) 

cat("Compiling forecasts", "\n")
newforecasts <- forecastall(forecast_date)
unlink("tmp/*")

cat("Rendering site", "\n")
rmarkdown::render_site()
