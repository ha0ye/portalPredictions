source('tools/model_functions.R')
source('tools/forecast_tools.R')
library(yaml)

model_metadata = yaml.load_file("data/model_metadata.yaml")
forecast_date = as.Date(model_metadata$forecast_date)


  # Get data

    all <- read.csv("data/rodent_all.csv")
    controls <- read.csv("data/rodent_controls.csv")
    model_metadata <- yaml.load_file("data/model_metadata.yaml")
    forecast_date <- as.Date(model_metadata$forecast_date)
    filename_suffix <- model_metadata$filename_suffix
    forecast_months <- model_metadata$forecast_months
    forecast_years <- model_metadata$forecast_years
    forecast_newmoons <- model_metadata$forecast_newmoons

    nm <- min(all$newmoonnumber):max(all$newmoonnumber)


    # grab the dates to go with count data

      moondat <- read.csv(text = RCurl::getURL(paste(
                            "https://raw.githubusercontent.com/",
                             "weecology/PortalData/master/Rodents/",
                             "moon_dates.csv", sep = "")),
                          stringsAsFactors = F)

      moondat$date <- as.Date(moondat$censusdate)
      nm_dates <- dplyr::filter(moondat, newmoonnumber %in% 
                                                nm) %>% 
                      dplyr::select(newmoonnumber, newmoondate)

    cdates <- as.Date(nm_dates$newmoondate)
    yr <- format(cdates, "%Y")

    foy <- round(
      as.numeric(format(cdates, "%j")) / 
      as.numeric(format(as.Date(paste(yr, "-12-31", sep = "")), "%j")), 3)


    Y <- rep(NA, length(nm))
   
    for(i in 1:length(nm)){

      ref <- which(all$newmoonnumber == nm[i])
      if(length(ref) > 0){
        Y[i] <- all$total[ref]
      }
    }

    data.frame(nm, Y, foy) 

nm <- nm[which(is.na(Y) == FALSE)]
foy <- foy[which(is.na(Y) == FALSE)]
Y <- Y[which(is.na(Y) == FALSE)]


  
