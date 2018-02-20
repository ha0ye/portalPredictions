
ft <- newmoons / max(newmoons)

b1 <- 0.0003
b2 <- 0.0019
b3 <- -0.0008
b4 <- 0.0005
b5 <- -0.00004
b6 <- 0.004

meanZ_a <- b1 * cos(2 * pi * frac_of_yr) + b2 * sin(2 * pi * frac_of_yr)
meanZ_d <- b3 * cos(2 * pi * frac_of_decade) + b4 * sin(2 * pi * frac_of_decade)
meanZ_y <- b5 * ft + b6 * ft^2
meanZ <- meanZ_a + meanZ_d + meanZ_y
tauZ <- 94

z <- rep(NA, length(Y) + 1)
z[1] <- 1.53
for(i in 1:length(Y)){
z[i + 1] <- z[i] + rnorm(1, meanZ[i], 1/sqrt(tauZ))
}

m4 <- apply(data.frame((exp(z) + 0.2), 0), 1, max)
m11 <- apply(data.frame((exp(z) - 0.14), 0), 1, max)
m14 <- apply(data.frame((exp(z) - 0.37), 0), 1, max)
m17 <- apply(data.frame((exp(z) + 0.31), 0), 1, max)


y4 <- rpois(length(z), m4)
y11 <- rpois(length(z), m11)
y14 <- rpois(length(z), m14)
y17 <- rpois(length(z), m17)

par(mfrow = c(2,1))
plot(y4, type = 'l', ylim = c(0, 60))
points(y11, type = 'l')
points(y14, type = 'l')
points(y17, type = 'l')


plot(t4, type = 'l', ylim = c(0,60))
points(t11, type = 'l')
points(t14, type = 'l')
points(t17, type = 'l')