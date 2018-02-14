source('tools/model_functions.R')
source('tools/forecast_tools.R')
library(yaml)

model_metadata <- yaml.load_file("data/model_metadata.yaml")
forecast_date <- as.Date(model_metadata$forecast_date)


plot_data <- portalr::abundance(level = "Plot", type = "Rodents", 
                                length = "all", incomplete = FALSE, 
                                time = "newmoon")
total <- apply(plot_data[ , 4:24], 1, sum)
plot_data <- data.frame(plot_data, total)

newmoons <- min(plot_data$newmoonnumber):max(plot_data$newmoonnumber)
n_newmoons <- length(newmoons)

moon_data <- get_moon_data()
moon_data$date <- as.Date(moon_data$newmoondate)
newmoon_dates <- dplyr::filter(moon_data, newmoonnumber %in% newmoons) %>% 
                      dplyr::select(newmoonnumber, newmoondate)

newmoon_date <- dates <- as.Date(newmoon_dates$newmoondate)
yr <- format(newmoon_date, "%Y")

newmoon_jday <- as.numeric(format(newmoon_date, "%j"))
nye_jday <- as.numeric(format(as.Date(paste(yr, "-12-31", sep = "")), "%j"))
frac_of_yr <- round(newmoon_jday / nye_jday, 3)

ref_1 <- paste(substr(yr, 1, 3), "0-01-01", sep = "")
ref_2 <- paste(substr(yr, 1, 3), "9-12-31", sep = "")

newmoon_days_from_ref_1 <- as.numeric(difftime(newmoon_date, ref_1))
ref_2_days_from_ref_1 <- as.numeric(difftime(ref_2, ref_1))
days_in_decade <- ref_2_days_from_ref_1 + 1

frac_of_decade <- round(newmoon_days_from_ref_1 / days_in_decade, 3)


   
specific_plot <- 4
specific_species <- "total"


Y <- rep(NA, n_newmoons)
Y_all <- matrix(NA, n_newmoons, length(4:25))

for(i in 1:n_newmoons){

  ref <- which(plot_data$newmoonnumber == newmoons[i] & 
               plot_data$plot == specific_plot)
  if(length(ref) > 0){
    Y[i] <- plot_data[ref, specific_species]
  }
}

for(i in 1:n_newmoons){

  ref <- which(plot_data$newmoonnumber == newmoons[i] & 
               plot_data$plot == specific_plot)
  if(length(ref) > 0){
    Y_all[i,] <- as.matrix(plot_data[ref, 4:25])
  }
}

windows(6.5, 8)

Y <- t4 + t11 + t14 + t17

par(mfrow = c(2, 1))
 
plot_title <- paste(specific_species, " in plot ", specific_plot, sep = "")
plot_title <- "total in plots 4, 11, 14, & 17 summed"
plot(newmoon_date, Y, type = 'l', las = 1, main = plot_title)
acf(Y, na.action = na.pass, lag.max = 100, las = 1, main = plot_title)




 
specific_plot <- 17
specific_species <- "total"


Y <- rep(NA, n_newmoons)
for(i in 1:n_newmoons){

  ref <- which(plot_data$newmoonnumber == newmoons[i] & 
               plot_data$plot == specific_plot)
  if(length(ref) > 0){
    Y[i] <- plot_data[ref, specific_species]
  }
}

t17 <- Y

ccf(t4, t11, na.action = na.pass, main = "total in plots 4 and 11", las = 1,
    lag.max = 25, ylim = c(-0.1, 1.0))

mean(c(cor(t4, t11, use = "complete.obs"),
cor(t4, t14, use = "complete.obs"),
cor(t4, t17, use = "complete.obs"),
cor(t11, t14, use = "complete.obs"),
cor(t11, t17, use = "complete.obs"),
cor(t14, t17, use = "complete.obs"))
)

tt <- data.frame(t4, t11, t14, t17)
mtt <- apply(tt, 1, mean)
vtt <- apply(tt, 1, var)

dvtt <- vtt/mtt - 1

hist(vtt/mtt, breaks = seq(0, 6, 0.1), las = 1, 
 main = "totals, plots 4, 11, 14, 17",
 xlab = "Variance/Mean")

plot(mtt, vtt, ylim = c(0, 50), xlim = c(0, 50), 
main = "totals, plots 4, 11, 14, 17", las = 1, ylab = "Variance", 
xlab = "Mean")
abline(0, 1)

plot(mtt, dvtt)

mm <- matrix(NA, nrow(tt), ncol(tt))
for(i in 1:nrow(tt)){
mm[i,] <- rpois(4, mtt[i])
}

rvtt <- apply(mm, 1, var)

hist(rvtt/mtt, breaks = seq(0, 6, 0.1))
plot(mtt, rvtt, ylim = c(0, 50), xlim = c(0, 50))
abline(0, 1)



plot(jitter(t4), jitter(t11), las = 1, xlab = "total in plot 4",
     ylab = "total in plot 11" )

plot(newmoon_date, Y, type = 'l',
     xlim = as.Date(c("1975-01-01", "1995-01-01")))
acf(Y[newmoons < "1995-01-01"], na.action = na.pass, lag.max = 50)

plot(newmoon_date, Y, type = 'l',
     xlim = as.Date(c("1995-01-01", "2018-01-01")))
acf(Y[newmoons >= "1995-01-01"], na.action = na.pass, lag.max = 50)

  
