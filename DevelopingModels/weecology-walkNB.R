#
# Weecology-WalkNB
#
# WalkNB is a model with a single mean density (on the log scale) that walks
#  through time. Sample variation is Negative Binomial
#

################
#
# it looks like the omega term explodes, which means the NB is just a 
#  Poisson, so this is on hold right now.
#
#  it's looking better now that i remembered what z0 means now
#    however it still explodes, which is why we should use precision
#  working with precision now, seems to be doing better (fingers crossed)
# it's not blowing up but i did get a fail with the library thing again in 
# mif2


library(pomp)
library(doParallel)

nm <- nm[which(is.na(Y) == FALSE)]
Y <- Y[which(is.na(Y) == FALSE)]

# model components 

  # process model

    process.sim <- "
      Z = Z + rnorm(0, 1 / sigmaZ);
    "

  # observation model

    count.rmeas <- "

      Y = rnbinom_mu(1 / omega, exp(Z));
    "

    count.dmeas <- "

      lik = dnbinom_mu(Y, 1 / omega, exp(Z), give_log);
    "

  # parameter transformations

    # to

      toES <- "
        Tomega = log(omega);
        TsigmaZ = log(sigmaZ);
      "

    # from

      fromES <- "
        Tomega = exp(omega);
        TsigmaZ = exp(sigmaZ);
      "

  mod <- pomp(data = data.frame(time = nm, Y = Y), 
               times = "time", t0 = 216,
               rprocess = discrete.time.sim(step.fun = Csnippet(process.sim),
                                            delta.t = 1), 
               rmeasure = Csnippet(count.rmeas),
               dmeasure = Csnippet(count.dmeas),
               paramnames = c("omega", "sigmaZ", 
                              "Z.0"),
               params = c(omega = 0.002, sigmaZ = 4.5,
                          Z.0 = 3.75),
               statenames = c("Z"),
               toEstimationScale = Csnippet(toES),	
               fromEstimationScale = Csnippet(fromES))

  plot(simulate(mod))
  pf <- pfilter(mod, Np = 1000)
  logLik(pf)


  Z0s <- seq(1, 5, 0.5)
  omegas <- seq(0.1, 5.1, 1)
  sigmaZs <- seq(10, 100, 10)
  ep <- expand.grid(Z0s, omegas, sigmaZs)
  colnames(ep) <- c("Z0", "omega", "sigmaZ")

  ncores <- 6
  doParallel::registerDoParallel(ncores)


lls <- foreach::foreach(i = 1:nrow(ep), .combine = rbind, 
                           .packages = c("pomp", "magrittr")) %dopar% {
    pfilter(mod, params = c(omega = ep$omega[i],
                             sigmaZ = ep$sigmaZ[i], 
                             Z.0 = ep$Z0[i]), 
             Np = 1000)$loglik
  }

stopImplicitCluster()

par(mfrow = c(3,1))
plot(ep$Z0, lls)
plot(ep$omega, lls)
plot(ep$sigmaZ, lls)






  Z0s <- seq(1, 5, 1)
  omegas <- seq(0.001, 0.011, 0.002)
  sigmaZs <- seq(1, 21, 5)
  guesses <- expand.grid(Z0s, omegas, sigmaZs)
  colnames(guesses) <- c("Z0", "omega", "sigmaZ")

  ncores <- 6
  doParallel::registerDoParallel(ncores)

  runs <- foreach::foreach(i = 1:nrow(guesses),  
                           .packages = c("pomp", "magrittr")) %dopar% {
                 mif2(mod, 
                      start = c("Z.0" = guesses[i, 1],
                                "omega" = guesses[i, 2],
                                "sigmaZ" = guesses[i, 3]),
                      Nmif = 200, Np = 1000, 
                      rw.sd = rw.sd(Z.0 = 0.02,
                                    omega = 0.01,
                                    sigmaZ = 0.08),
                      cooling.fraction.50 = 0.95, 
                      cooling.type = "hyperbolic",
                      transform = TRUE)  %>% 
                 continue(Nmif = 200, cooling.fraction = 0.8) %>%
                 continue(Nmif = 200, cooling.fraction = 0.6) %>%
                 continue(Nmif = 200, cooling.fraction = 0.4) %>%
                 continue(Nmif = 200, cooling.fraction= 0.2) -> m1
                 m1
  }





windows(12,12)
par(mfrow = c(3, 1))
plot(1, 1, type = 'n', ylim = (c(0, 5)), xlim = c(0, 1000), 
     ylab = "Z.0", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"Z.0"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}

plot(1, 1, type = 'n', ylim = (c(0,0.01)), xlim = c(0, 1000), 
     ylab = "omega", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"omega"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}

plot(1, 1, type = 'n', ylim = (c(3, 6)), xlim = c(0, 1000), 
     ylab = "sigmaZ", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"sigmaZ"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}



