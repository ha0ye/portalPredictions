#
# for the ESA abstract
#
#
#  things to do: 
#    set up an initialization function so you can more properly model things with hyper parameters!



#
library(pomp)
library(doParallel)


# model components 

  # process model

    process.sim <- "
      double meanZ_a = b1 * cos(M_2PI * foy) + b2 * sin(M_2PI * foy);
      double meanZ_d = b3 * cos(M_2PI * fod) + b4 * sin(M_2PI * fod);
      double meanZ_y = b5 * ft + b6 * ft * ft;
      double meanZ = meanZ_y + meanZ_a + meanZ_d;
      double sigmaZ = 1 / sqrt(tauZ);
      Z = Z + rnorm(meanZ, sigmaZ);
    "

  # observation model

    count.rmeas <- "

      double m4 = fmax2(exp(Z) + 0.20, 0);
      double m11 = fmax2(exp(Z) - 0.14, 0);
      double m14 = fmax2(exp(Z) - 0.37, 0);
      double m17 = fmax2(exp(Z) + 0.31, 0);
      int p4 = rpois(m4);
      int p11 = rpois(m11);
      int p14 = rpois(m14);
      int p17 = rpois(m17);
      t4 = rbinom(p4, rho);
      t11 = rbinom(p11, rho);
      t14 = rbinom(p14, rho);
      t17 = rbinom(p17, rho);
    "

    count.dmeas <- "

      double m4 = fmax2(exp(Z) + 0.20, 0);
      double m11 = fmax2(exp(Z) - 0.14, 0);
      double m14 = fmax2(exp(Z) - 0.37, 0);
      double m17 = fmax2(exp(Z) + 0.31, 0);
      double p4 = rpois(m4);
      double p11 = rpois(m11);
      double p14 = rpois(m14);
      double p17 = rpois(m17);
      lik = dbinom(t4, p4, rho, give_log) + 
            dbinom(t11, p11, rho, give_log) + 
            dbinom(t14, p14, rho, give_log) + 
            dbinom(t17, p17, rho, give_log);
    "

  # parameter transformations

    # to

      toES <- "
        TtauZ = log(tauZ);
      "

    # from

      fromES <- "
        TtauZ = exp(tauZ);
      "

newmoons_nar <- newmoons[-which(is.na(t4) == TRUE)]
frac_of_yr_nar <- frac_of_yr[-which(is.na(t4) == TRUE)]
frac_of_decade_nar <- frac_of_decade[-which(is.na(t4) == TRUE)]
t4_nar <- t4[-which(is.na(t4) == TRUE)]
t11_nar <- t11[-which(is.na(t4) == TRUE)]
t14_nar <- t14[-which(is.na(t4) == TRUE)]
t17_nar <- t17[-which(is.na(t4) == TRUE)]

    mod <- pomp(data = data.frame(time = newmoons_nar, t4 = t4_nar,
                                  t11 = t11_nar, t14 = t14_nar, t17 = t17_nar), 
               times = "time", t0 = 0,
              covar = data.frame(time = c(0, newmoons_nar), 
                        ft = c(0, newmoons_nar) / max(newmoons_nar), 
                        foy = c(frac_of_yr_nar[1] - 0.07, frac_of_yr_nar), 
                        fod = c(frac_of_decade_nar[1] - 0.07, frac_of_decade_nar)), 
               tcovar = "time",
               rprocess = discrete.time.sim(step.fun = Csnippet(process.sim),
                                            delta.t = 1), 
               rmeasure = Csnippet(count.rmeas),
               dmeasure = Csnippet(count.dmeas),
               paramnames = c("b1", "b2", "b3", "b4", "b5", "b6", "tauZ", 
                              "Z.0",
                              "rho"),
               params = c(b1 = 0.00056, b2 = 0.0018, 
                          b3 = -0.0014, b4 = 0.000057, 
                          b5 = 0.0004, b6 = -0.00056, 
                          tauZ = 97.8,
                          Z.0 = 1.66,
                          rho = 0.781),
               statenames = c("Z"),
               toEstimationScale = Csnippet(toES),	
               fromEstimationScale = Csnippet(fromES))


  pf <- pfilter(mod, Np = 1000)
  logLik(pf)

ff<-(simulate(mod))
plot(ff)
ccf(t(ff@data)[,1],t(ff@data)[,2])
mean(apply(ff@data, 2, var) / apply(ff@data, 2, mean))



b3 = 0.003, b4 = -0.005
b5 = 0.002, b6 = 0.005

  Z0s <- runif(5000, 0.5, 2)
  b1s <- runif(5000, -0.005, 0.005)
  b2s <- runif(5000, -0.005, 0.005)
  b3s <- runif(5000, -0.005, 0.005)
  b4s <- runif(5000, -0.005, 0.005)
  b5s <- runif(5000, -0.005, 0.005)
  b6s <- runif(5000, -0.005, 0.005)
  tauZs <- runif(5000, 10, 1000)
  rhos <- runif(5000, 0.5, 1)

  ep <- data.frame(Z0s, b1s, b2s, b3s, b4s, b5s, b6s, tauZs, rhos)
  colnames(ep) <- c("Z0", "b1", "b2", "b3", "b4", "b5", "b6", "tauZ", "rho")

  ncores <- 6
  doParallel::registerDoParallel(ncores)

  lls <- foreach::foreach(i = 1:nrow(ep), .combine = rbind, .errorhandling = "pass", 
                           .packages = c("pomp", "magrittr")) %dopar% {
     pfilter(mod, params = c(tauZ = ep$tauZ[i], b1 = ep$b1[i], b2 = ep$b2[i], 
                              b3 = ep$b3[i], b4 = ep$b4[i],  b5 = ep$b5[i], b6 = ep$b6[i], 
                             Z.0 = ep$Z0[i], rho = ep$rho[i]), 
             Np = 1000)$loglik
  }

stopImplicitCluster()


par(mfrow = c(5,2), mar = c(5,4,1,1))
plot(ep$Z0, lls, col = rgb(0.5, 0.5, 0.5, 0.3))
plot(ep$tauZ, lls, col = rgb(0.5, 0.5, 0.5, 0.3))
plot(ep$b1, lls, col = rgb(0.5, 0.5, 0.5, 0.3))
plot(ep$b2, lls, col = rgb(0.5, 0.5, 0.5, 0.3))
plot(ep$b3, lls, col = rgb(0.5, 0.5, 0.5, 0.3))
plot(ep$b4, lls, col = rgb(0.5, 0.5, 0.5, 0.3))
plot(ep$b5, lls, col = rgb(0.5, 0.5, 0.5, 0.3))
plot(ep$b6, lls, col = rgb(0.5, 0.5, 0.5, 0.3))
plot(ep$rho, lls, col = rgb(0.5, 0.5, 0.5, 0.3))


z0 .5 - 2, 1.5
b1 -0.002
b2 0.0606
b3 -0.0019
b4 0.0039
b5 
b6 0.006


  Z0s <- runif(18, 1, 5)
  b1s <- runif(18, -3, 3)
  b2s <- runif(18, -3, 3)
  b3s <- runif(18, -3, 3)
  b4s <- runif(18, -3, 3)
  b5s <- runif(18, -3, 3)
  b6s <- runif(18, -3, 3)
  tauZs <- runif(12, exp(-5), exp(-0.5))
  guesses <- data.frame(Z0s, tauZs, b1s, b2s, b3s, b4s, b5s, b6s)
  colnames(guesses) <- c("Z0", "tauZ", "b1", "b2", "b3", "b4", "b5", "b6")


  ncores <- 6
  doParallel::registerDoParallel(ncores)

  runs <- foreach::foreach(i = 1:nrow(guesses),  
                           .packages = c("pomp", "magrittr")) %dopar% {
                 mif2(mod, 
                      start = c("Z.0" = guesses[i, 1],
                                "sigmaZ" = guesses[i, 2],
                                "b1" = guesses[i, 3],
                                "b2" = guesses[i, 4],
                                "b3" = guesses[i, 5],
                                "b4" = guesses[i, 6],
                                "b5" = guesses[i, 7],
                                "b6" = guesses[i, 8]),
                      Nmif = 100, Np = 1000, 
                      rw.sd = rw.sd(Z.0 = 0.01,
                                    sigmaZ = 0.002,
                                    b1 = 0.001,
                                    b2 = 0.001,
                                    b3 = 0.001,
                                    b4 = 0.001,
                                    b5 = 0.001,
                                    b6 = 0.001),
                      cooling.fraction.50 = 0.25, 
                      cooling.type = "geometric",
                      transform = TRUE)  %>% 
                 continue(Nmif = 100, cooling.fraction = 0.5) %>%
                 continue(Nmif = 100, cooling.fraction = 0.75) %>%
                 continue(Nmif = 100, cooling.fraction = 0.95) -> m1
                 m1
  }

stopImplicitCluster()






par(mfrow = c(4, 2))
plot(1, 1, type = 'n', ylim = (c(-2, 6)), xlim = c(0, 400), 
     ylab = "Z.0", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"Z.0"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}
plot(1, 1, type = 'n', ylim = (c(0, 1)), xlim = c(0, 400), 
     ylab = "sigmaZ", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"sigmaZ"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}
plot(1, 1, type = 'n', ylim = (c(-1, 1)), xlim = c(0, 400), 
     ylab = "b1", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"b1"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}
plot(1, 1, type = 'n', ylim = (c(-1, 1)), xlim = c(0, 400), 
     ylab = "b2", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"b2"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}

plot(1, 1, type = 'n', ylim = (c(-1, 1)), xlim = c(0, 400), 
     ylab = "b3", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"b3"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}
plot(1, 1, type = 'n', ylim = (c(-1, 1)), xlim = c(0, 400), 
     ylab = "b4", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"b4"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}

plot(1, 1, type = 'n', ylim = (c(-1, 1)), xlim = c(0, 400), 
     ylab = "b5", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"b5"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}
plot(1, 1, type = 'n', ylim = (c(-1, 1)), xlim = c(0, 400), 
     ylab = "b6", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"b6"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}


















Z0_mle <- rep(NA, 18)
for(i in 1:18){
Z0_mle[i] <- conv.rec(runs[[i]])[400,"Z.0"]
}

sigmaZ_mle <- rep(NA, 18)
for(i in 1:18){
sigmaZ_mle[i] <- conv.rec(runs[[i]])[400,"sigmaZ"]
}

b1_mle <- rep(NA, 18)
for(i in 1:18){
b1_mle[i] <- conv.rec(runs[[i]])[400,"b1"]
}

b2_mle <- rep(NA, 18)
for(i in 1:18){
b2_mle[i] <- conv.rec(runs[[i]])[400,"b2"]
}

par(mfrow = c(4,1))
hist(Z0_mle, breaks = seq(3, 6, 0.1))
hist(sigmaZ_mle, breaks = seq(0, 0.5, 0.01))
hist(b1_mle, breaks = seq(-0.1, 0.1, 0.01))
hist(b2_mle, breaks = seq(0, 0.2, 0.01))