
library(portalr)



?portalr::get_rodent_data


x <- portalr::abundance(level="Plot",type="Rodents",length="all", 
                        incomplete = FALSE)

xx<-x[x$treatment == "control",]
table(xx$plot)





plots
1 9
2 8 12 22
4 11 14 17




source('tools/model_functions.R')

moon_data <- get_moon_data()

rodent_data <- portalr::get_rodent_data(path = "repo", level = "Plot", 
                                        length = "all", time = "newmoon")

total <- apply(rodent_data[, 4:24], 1, sum)
rodent_data <- data.frame(rodent_data, total)

moons <- moon_data$newmoonnumber
n_moons <- length(moons)

plots <- unique(rodent_data$plot)
n_plots <- length(plots)

moons_plots <- expand.grid(moon = moons, plot = plots)
moonplot <- paste(moons_plots$moon, moons_plots$plot, sep = "-")

rd_moonplot <- paste(rodent_data$newmoonnumber, rodent_data$plot, sep = "-")

which_in <- which(moonplot %in% rd_moonplot)
which_not_in <- (1:length(moonplot))[-which_in]

to_add <- matrix(NA, nrow = length(which_not_in), ncol = ncol(rodent_data))
to_add[ , 1:2] <- as.matrix(moons_plots[which_not_in, ])
colnames(to_add) <- colnames(rodent_data)

rodent_data_expanded <- rbind(rodent_data, to_add)


nobs <- nrow(rodent_data_expanded)
census_date <- rep(NA, nobs)
newmoon_date <- rep(NA, nobs)
for(i in 1:nobs){
  which_obs <- which(moon_data$newmoonnumber == 
                       rodent_data_expanded$newmoonnumber[i])
  census_date[i] <- as.character(moon_data$censusdate[which_obs])
  newmoon_date[i] <- as.character(moon_data$newmoondate[which_obs])
}
census_date <- as.Date(census_date)
newmoon_date <- as.Date(newmoon_date)

rodent_data_expanded <- data.frame(rodent_data_expanded, 
                                   censusdate = census_date, 
                                   newmoondate = newmoon_date)

censuses_to_use <- which(rodent_data_expanded$newmoonnumber >= 4)
rodent_data_to_use <- rodent_data_expanded[censuses_to_use, ]

reorder <- order(rodent_data_to_use$plot, rodent_data_to_use$newmoonnumber)
rodent_data_to_use <- rodent_data_to_use[reorder, ]


# exclosure = k rats removed (DM, DO, DS)
# removal = all rodents removed

unique(rodent_data_to_use$plot[which(
           rodent_data_to_use$treatment == "spectabs")])


specific_plot <- 11
specific_sp <- "DM"
tp <- rodent_data_to_use[rodent_data_to_use$plot == specific_plot, 
                         specific_sp]

par(mfrow = c(2,2))
plot_title <- paste(specific_plot, specific_sp, sep = " ")
plot(tp, type = 'l', main = plot_title)
acf(tp, na.action = na.pass, lag.max = 200, main = plot_title)
plot_title <- paste(specific_plot, specific_sp, "second half", sep = " ")
plot(tp[250:500], type = 'l', main = plot_title)
acf(tp[250:500], na.action = na.pass, lag.max = 100, main = plot_title)

