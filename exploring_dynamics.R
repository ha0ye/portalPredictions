

    fod <- round(
 ( as.numeric(difftime(cdates, "1990-01-01")) %% (10*365.3) ) / (10*365.3), 3)

fod2 <- c(0.492, fod)
foy2 <- c(0.923, foy)
z <- rep(NA, length(foy2))
y <- rep(NA, length(foy2))

z0 <- 4
z[1] <- z0
b1 <- 0.35
b2 <- 0.05
bd1 <- -0.015
bd2 <- -0.0
r <- b1 * cos(2 * pi * foy2) + b2 * sin(2 * pi * foy2)
rd <- bd1 * cos(2 * pi * fod2) + bd2 * sin(2 * pi * fod2)
sig <- 0.03
sig2 <- 0.0

for(i in 2:length(z)){
  z[i] <- z[i - 1] + rnorm(1, 0, sig) + r[i - 1] + rd[i - 1]
  y[i] <- rpois(1, exp(z[i]+ rnorm(1, 0, sig2)))
}

plot(y, type = 'l')
points(y, col = 3)