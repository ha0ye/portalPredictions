#
# Weecology-SeasonalWalk
#
# SeasonalWalk is a model with a single mean density (on the log scale) that 
#  walks and is seasonally dynamic through time. Sample variation is Poisson
#
#
# very much in the thick of things here. trying to figure out how best to incorporate 
#  the seasonal dynamics to match the data
#
#
library(pomp)
library(doParallel)


# model components 

  # process model

    process.sim <- "
      Z = Z + Z * rnorm(b1 * cos(M_2PI * foy) + b2 * sin(M_2PI * foy), sigmaZ);
    "

  # observation model

    count.rmeas <- "

      Y = rpois(exp(Z));
    "

    count.dmeas <- "

      lik = dpois(Y, exp(Z), give_log);
    "

  # parameter transformations

    # to

      toES <- "
        TsigmaZ = log(sigmaZ);
      "

    # from

      fromES <- "
        TsigmaZ = exp(sigmaZ);
      "

  mod <- pomp(data = data.frame(time = nm, Y = Y), 
               times = "time", t0 = 216,
              covar = data.frame(time = c(216, nm), 
                                 foy = c(foy[1] - 0.07, foy)), 
               tcovar = "time",
               rprocess = discrete.time.sim(step.fun = Csnippet(process.sim),
                                            delta.t = 1), 
               rmeasure = Csnippet(count.rmeas),
               dmeasure = Csnippet(count.dmeas),
               paramnames = c("b1", "b2", "sigmaZ", 
                              "Z.0"),
               params = c(b1 = 0.005, b2 = 0.035, sigmaZ = 0.025,
                          Z.0 = 3.79),
               statenames = c("Z"),
               toEstimationScale = Csnippet(toES),	
               fromEstimationScale = Csnippet(fromES))


  pf <- pfilter(mod, Np = 1000)
  logLik(pf)


  Z0s <- runif(100, 1, 5)
  b1s <- runif(100, -2, 2)
  b2s <- runif(100, -2, 2)
  sigmaZs <- runif(100, exp(-3), exp(0))
  ep <- data.frame(Z0s, b1s, b2s, sigmaZs)
  colnames(ep) <- c("Z0", "b1", "b2", "sigmaZ")

  ncores <- 6
  doParallel::registerDoParallel(ncores)

  lls <- foreach::foreach(i = 1:nrow(ep), .combine = rbind, 
                           .packages = c("pomp", "magrittr")) %dopar% {
     pfilter(mod, params = c(sigmaZ = ep$sigmaZ[i], b1 = ep$b1[i], b2 = ep$b2[i], 
                             Z.0 = ep$Z0[i]), 
             Np = 1000)$loglik
  }

stopImplicitCluster()
par(mfrow = c(2,1))
plot(ep$Z0, lls)
plot(ep$sigmaZ, lls)

par(mfrow = c(4,1))
plot(ep$Z0, lls)
plot(ep$b1, lls)
plot(ep$b2, lls)
plot(ep$sigmaZ, lls)


  Z0s <- runif(18, 1, 5)
  b1s <- runif(18, -1, 1)
  b2s <- runif(18, -1, 1)
  sigmaZs <- runif(18, 0.01, 1)
  guesses <- data.frame(Z0s, b1s, b2s, sigmaZs)
  colnames(guesses) <- c("Z0", "b1", "b2", "sigmaZ")

  ncores <- 6
  doParallel::registerDoParallel(ncores)

  runs <- foreach::foreach(i = 1:nrow(guesses),  
                           .packages = c("pomp", "magrittr")) %dopar% {
                 mif2(mod, 
                      start = c("Z.0" = guesses[i, 1],
                                "sigmaZ" = guesses[i, 4],
                                "b1" = guesses[i, 2],
                                "b2" = guesses[i, 3]),
                      Nmif = 100, Np = 1000, 
                      rw.sd = rw.sd(Z.0 = 0.01,
                                    sigmaZ = 0.002,
                                    b1 = 0.001,
                                    b2 = 0.001),
                      cooling.fraction.50 = 0.25, 
                      cooling.type = "geometric",
                      transform = TRUE) %>% 
                 continue(Nmif = 100, cooling.fraction = 0.5) %>%
                 continue(Nmif = 100, cooling.fraction = 0.75) %>%
                 continue(Nmif = 100, cooling.fraction = 0.95)  -> m1
                 m1
  }

stopImplicitCluster()



par(mfrow = c(4, 1))
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