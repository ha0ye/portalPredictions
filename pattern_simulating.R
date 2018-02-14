b0 <- 0.2
b1 <- -0.4
b2 <- 0.0 
sig <- 0.05
z1 <- b0 + b1 * cos(2 * pi * frac_of_yr) + b2 * sin(2 * pi * frac_of_yr)
er1 <- rep(NA, length(z1))
val1 <- 0
for(i in 1:length(z1)){
  val1 <- val1 + rnorm(1, 0, sig)
  er1[i] <- val1
}


b0 <- 0.2
b1 <- -0.
b2 <- 0.1 
sig <- 0.2
z2 <- b0 + b1 * cos(2 * pi * frac_of_yr) + b2 * sin(2 * pi * frac_of_yr)
er2 <- rep(NA, length(z2))
val2 <- 0
for(i in 1:length(z2)){
  val2 <- val2 + rnorm(1, 0, sig)
  er2[i] <- val2
}

zt <- z1 + z2
et <- er1 + er2

zzt <- rpois(length(zt), exp(zt + et))
par(mfrow = c(2,2))
plot(newmoon_date, zzt, type = 'l')
acf(zzt, lag.max = 100)
plot(frac_of_yr, exp(z1))
points(frac_of_yr, exp(z2))
plot(newmoon_date, er2, type = 'l')
points(newmoon_date, er1, type = 'l')








