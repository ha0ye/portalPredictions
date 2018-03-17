source("tools/forecast_tools.R")

model_metadata <- yaml::yaml.load_file("data/model_metadata.yaml")
forecast_date <- as.Date(model_metadata$forecast_date)

if(!(forecast_date == Sys.Date())){ 
  stop('Data not updated') 
}

# Run all models 
#  this will need to be generalized for models written in other languages
cat("Running models", "\n")
dir.create("tmp")
sapply(list.files("models", full.names = TRUE), source) 

# Collect all forecast results and save to predictions directory
cat("Compiling forecasts", "\n")
newforecasts <- forecastall(forecast_date)
unlink("tmp/*")

# Update Website
rmarkdown::render_site()
