context("checks that new forecasts are being added correctly")
model_metadata <- yaml::yaml.load_file("../../data/model_metadata.yaml")

forecast_date <- as.Date(model_metadata$forecast_date)
file_suffix <- model_metadata$filename_suffix
forecast_months <- model_metadata$rodent_forecast_months
forecast_years <- model_metadata$rodent_forecast_years
forecast_newmoons <- model_metadata$rodent_forecast_newmoons

fcast_fname <- file.path('../../predictions', 
                 paste(forecast_date, file_suffix, ".csv", sep=""))
forecasts <- read.csv(fcast_fname, na.strings = "")

aic_fname <- file.path('../../predictions', 
                 paste(forecast_date, file_suffix, "_model_aic.csv", sep=""))
forecastaics <- read.csv(aic_fname, na.strings = "")

forecastnames <- c("date", "forecastmonth", "forecastyear", "newmoonnumber",
                   "currency", "model", "level", "species", "estimate",  
                   "LowerPI", "UpperPI", "fit_start_newmoon", 
                   "fit_end_newmoon", "initial_newmoon")

aicnames <- c("date", "currency", "model", "level", "species", "aic", 
              "fit_start_newmoon", "fit_end_newmoon", "initial_newmoon")

valid_currencies <- c("abundance", "biomass","energy")
valid_levels <- paste("Plot", 1:24, " ", sep = "")
valid_levels <- c("All", "Controls", "FullExclosure", "KratExclosure", 
                  valid_levels)
valid_species <- c("total", "BA", "DM", "DO", "DS", "OL", "OT", "PB", "PE", 
                   "PF", "PH", "PL", "PM", "PP", "RF", "RM", "RO", "SF", "SH",
                   "SO", "NA")

testthat::test_that("column names correct", {
  testthat::expect_true(all(colnames(forecasts) == forecastnames))
  testthat::expect_true(all(colnames(forecastaics) == aicnames))
})

testthat::test_that("dates are valid", {
  testthat::expect_false(any(is.na(as.Date(forecasts$date, "%Y-%m-%d"))))
  testthat::expect_false(any(is.na(as.Date(forecastaics$date, "%Y-%m-%d"))))
})

testthat::test_that("years and months are valid", {
  testthat::expect_true(
    all(forecasts$month %in% model_metadata$rodent_forecast_months))
  testthat::expect_true(
    all(forecasts$year %in% model_metadata$rodent_forecast_years))
})

testthat::test_that("newmoons are valid", {
  testthat::expect_true(
    all(forecasts$newmoonnumber %in% model_metadata$rodent_forecast_newmoons))
})

testthat::test_that("currencies are valid", {
  testthat::expect_true(all(forecasts$currency %in% valid_currencies))
  testthat::expect_true(all(forecastaics$currency %in% valid_currencies))
})

testthat::test_that("level is valid", {
  testthat::expect_true(all(forecasts$level %in% valid_levels))
  testthat::expect_true(all(forecastaics$level %in% valid_levels))
})

testthat::test_that("species is valid", {
  testthat::expect_true(all(forecasts$species %in% valid_species))
  testthat::expect_true(all(forecastaics$species %in% valid_species))
})

testthat::test_that("no NAs in estimates, PIs, or AICs", {
  testthat::expect_false(any(is.na(forecasts$estimate))) 
  testthat::expect_false(any(is.na(forecasts$LowerPI))) 
  testthat::expect_false(any(is.na(forecasts$UpperPI)))
  testthat::expect_false(any(is.na(forecastaics$aic)))
})

testthat::test_that("newmoon columns are integers only", {
  testthat::expect_true(is.integer(forecasts$fit_start_newmoon)) 
  testthat::expect_true(is.integer(forecasts$fit_end_newmoon)) 
  testthat::expect_true(is.integer(forecasts$initial_newmoon)) 
  testthat::expect_false(any(is.na(forecasts$fit_start_newmoon))) 
  testthat::expect_false(any(is.na(forecasts$fit_end_newmoon))) 
  testthat::expect_false(any(is.na(forecasts$initial_newmoon))) 
})
  
testthat::test_that("no forecasts duplicated", {
  testthat::expect_true(sum(duplicated(forecasts[,1:8])) == 0)
})