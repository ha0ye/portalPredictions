# i'm currently working on tidying the forecast_tools script

# here are libraries that i've cut out the calls to. need to use ::

library(tidyverse)
library(lubridate)
library(zoo)
library(ggplot2)
library(rmarkdown)
library(RCurl)


# notes cut from functions

  # make_ensemble

  #Mean is the weighted mean of all model means.
  #Variance is the weighted mean of all model variances + the variances of the weighted mean 
  #using the unbiased estimate of sample variance. See https://github.com/weecology/portalPredictions/pull/65
  #We only store the prediction interval for models, so backcalculate individual model variance
  #assuming the same CI_level throughout. 
  #Assert that the summed weight of all the model ensembles is 1, as that's what the above variance estimates assume.
  #Rounded to account for precision errors. Summed weights can also be NA if there are not weights availble for that ensemble. 
