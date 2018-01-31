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

      Y = rnbinom_mu(omega, exp(Z));
    "

    count.dmeas <- "

      lik = dnbinom_mu(Y, omega, exp(Z), give_log);
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
               params = c(omega = 1.0, sigmaZ = 0.24,
                          Z.0 = 3.5),
               statenames = c("Z"),
               toEstimationScale = Csnippet(toES),	
               fromEstimationScale = Csnippet(fromES))


  pf <- pfilter(mod, Np = 1000)
  logLik(pf)


  Z0s <- seq(1, 5, 0.5)
  omegas <- exp(seq(1, 3.5, 0.5))
  sigmaZs <- exp(seq(-3, 0, 0.5))
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

par(mfrow = c(3,1))
plot(ep$Z0, lls)
plot(ep$omega, lls)
plot(ep$sigmaZ, lls)






  Z0s <- seq(2, 5, 1)
  omegas <- seq(1.1, 5.1, 2)
  sigmaZs <- exp(seq(-2, 0.5, 0.5))
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
                      Nmif = 100, Np = 1000, 
                      rw.sd = rw.sd(Z.0 = 0.02,
                                    omega = 0.01,
                                    sigmaZ = 0.08),
                      cooling.fraction.50 = 0.95, 
                      cooling.type = "hyperbolic",
                      transform = TRUE)  -> m1
                 m1
  }


%>% 
                 continue(Nmif = 100, cooling.fraction = 0.8) %>%
                 continue(Nmif = 100, cooling.fraction = 0.6) %>%
                 continue(Nmif = 100, cooling.fraction = 0.4) %>%
                 continue(Nmif = 100, cooling.fraction= 0.2)


windows(12,12)
par(mfrow = c(3, 1))
plot(1, 1, type = 'n', ylim = exp(c(-2, 5)), xlim = c(0, 500), 
     ylab = "Z.0", xlab = "iteration")
for(i in 1:nrow(guesses)){
points(exp(conv.rec(runs[[i]])[,"Z.0"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}

plot(1, 1, type = 'n', ylim = (c(0,2000)), xlim = c(0, 500), 
     ylab = "omega", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"omega"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}

plot(1, 1, type = 'n', ylim = (c(0, 0.25)), xlim = c(0, 500), 
     ylab = "sigmaZ", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"sigmaZ"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}



