#
# Weecology-FixedNB
#
# FixedNB is a model with a single mean density (on the log scale) that does
#  not evolve over time. Sample variation is Negative Binomial
#

library(pomp)
library(doParallel)

nm <- nm[which(is.na(Y) == FALSE)]
Y <- Y[which(is.na(Y) == FALSE)]

# model components 

  # process model

    process.sim <- "
      Z = Z;
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
      "

    # from

      fromES <- "
        Tomega = exp(omega);
      "

  mod <- pomp(data = data.frame(time = nm, Y = Y), 
               times = "time", t0 = 216,
               rprocess = discrete.time.sim(step.fun = Csnippet(process.sim),
                                            delta.t = 1), 
               rmeasure = Csnippet(count.rmeas),
               dmeasure = Csnippet(count.dmeas),
               paramnames = c("omega",
                              "Z.0"),
               params = c(omega = 1.0, 
                          Z.0 = log(mean(na.omit(Y)))),
               statenames = c("Z"),
               toEstimationScale = Csnippet(toES),	
               fromEstimationScale = Csnippet(fromES))


  pf <- pfilter(mod, Np = 1000)
  logLik(pf)


  Z0s <- seq(4, 6, 0.25)
  omegas <- exp(seq(-0, 2, 0.25))
  ep <- expand.grid(Z0s, omegas)
  colnames(ep) <- c("Z0", "omega")

  ncores <- 6
  doParallel::registerDoParallel(ncores)

  lls <- foreach::foreach(i = 1:nrow(ep), .combine = rbind, 
                           .packages = c("pomp", "magrittr")) %dopar% {
     pfilter(mod, params = c(omega = ep$omega[i],
                             Z.0 = ep$Z0[i]), 
             Np = 1000)$loglik
  }

par(mfrow = c(2,1))
plot(ep$Z0, lls)
plot(ep$omega, lls)
plot(ep$omega[ep$Z0==5.25], lls[ep$Z0==5.25])



  Z0s <- seq(4, 6, 1)
  omegas <- seq(1.1, 5.1, 2)
  guesses <- expand.grid(Z0s, omegas)
  colnames(guesses) <- c("Z0", "omega")

  ncores <- 6
  doParallel::registerDoParallel(ncores)

  runs <- foreach::foreach(i = 1:nrow(guesses),  
                           .packages = c("pomp", "magrittr")) %dopar% {
                 mif2(mod, 
                      start = c("Z.0" = guesses[i, 1],
                                "omega" = guesses[i, 2]),
                      Nmif = 100, Np = 1000, 
                      rw.sd = rw.sd(Z.0 = 0.1,
                                    omega = 0.05),
                      cooling.fraction.50 = 0.95, 
                      cooling.type = "hyperbolic",
                      transform = TRUE) %>% 
                 continue(Nmif = 100, cooling.fraction = 0.8) %>%
                 continue(Nmif = 100, cooling.fraction = 0.6) %>%
                 continue(Nmif = 100, cooling.fraction = 0.4) %>%
                 continue(Nmif = 100, cooling.fraction= 0.2) -> m1
                 m1
  }


windows(12,12)
par(mfrow = c(2, 1))
plot(1, 1, type = 'n', ylim = exp(c(-2, 6)), xlim = c(0, 500), 
     ylab = "Z.0", xlab = "iteration")
for(i in 1:nrow(guesses)){
points(exp(conv.rec(runs[[i]])[,"Z.0"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}

plot(1, 1, type = 'n', ylim = exp(c(0,4)), xlim = c(0, 500), 
     ylab = "omega", xlab = "iteration")
for(i in 1:nrow(guesses)){
points(exp(conv.rec(runs[[i]])[,"omega"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}




