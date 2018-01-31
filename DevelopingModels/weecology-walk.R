#
# Weecology-Walk
#
# Walk is a model with a single mean density (on the log scale) that walks
#  through time. Sample variation is Poisson
#

library(pomp)
library(doParallel)

nm <- nm[which(is.na(Y) == FALSE)]
Y <- Y[which(is.na(Y) == FALSE)]

# model components 

  # process model

    process.sim <- "
      Z = Z + rnorm(0, sigmaZ);
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
               rprocess = discrete.time.sim(step.fun = Csnippet(process.sim),
                                            delta.t = 1), 
               rmeasure = Csnippet(count.rmeas),
               dmeasure = Csnippet(count.dmeas),
               paramnames = c("sigmaZ", 
                              "Z.0"),
               params = c(sigmaZ = 0.24,
                          Z.0 = 3.5),
               statenames = c("Z"),
               toEstimationScale = Csnippet(toES),	
               fromEstimationScale = Csnippet(fromES))


  pf <- pfilter(mod, Np = 1000)
  logLik(pf)


  Z0s <- seq(1, 5, 0.5)
  sigmaZs <- exp(seq(-3, 0, 0.5))
  ep <- expand.grid(Z0s, sigmaZs)
  colnames(ep) <- c("Z0", "sigmaZ")

  ncores <- 6
  doParallel::registerDoParallel(ncores)

  lls <- foreach::foreach(i = 1:nrow(ep), .combine = rbind, 
                           .packages = c("pomp", "magrittr")) %dopar% {
     pfilter(mod, params = c(sigmaZ = ep$sigmaZ[i], 
                             Z.0 = ep$Z0[i]), 
             Np = 1000)$loglik
  }

stopImplicitCluster()
par(mfrow = c(2,1))
plot(ep$Z0, lls)
plot(ep$sigmaZ, lls)






  Z0s <- seq(2, 5, 0.5)
  sigmaZs <- exp(seq(-2, 1, 0.5))
  guesses <- expand.grid(Z0s, sigmaZs)
  colnames(guesses) <- c("Z0", "sigmaZ")

  ncores <- 6
  doParallel::registerDoParallel(ncores)

  runs <- foreach::foreach(i = 1:nrow(guesses),  
                           .packages = c("pomp", "magrittr")) %dopar% {
                 mif2(mod, 
                      start = c("Z.0" = guesses[i, 1],
                                "sigmaZ" = guesses[i, 2]),
                      Nmif = 200, Np = 1000, 
                      rw.sd = rw.sd(Z.0 = 0.02,
                                    sigmaZ = 0.08),
                      cooling.fraction.50 = 0.95, 
                      cooling.type = "hyperbolic",
                      transform = TRUE) %>% 
                 continue(Nmif = 200, cooling.fraction = 0.8) %>%
                 continue(Nmif = 200, cooling.fraction = 0.6) %>%
                 continue(Nmif = 200, cooling.fraction = 0.4) %>%
                 continue(Nmif = 200, cooling.fraction= 0.2)  -> m1
                 m1
  }

stopImplicitCluster()




windows(12,12)
par(mfrow = c(2, 1))
plot(1, 1, type = 'n', ylim = exp(c(-2, 5)), xlim = c(0, 1000), 
     ylab = "Z.0", xlab = "iteration")
for(i in 1:nrow(guesses)){
points(exp(conv.rec(runs[[i]])[,"Z.0"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}
plot(1, 1, type = 'n', ylim = (c(0.1, 0.3)), xlim = c(0, 1000), 
     ylab = "sigmaZ", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"sigmaZ"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}



