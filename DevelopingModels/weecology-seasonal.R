#
# Weecology-Seasonal
#
# Seasonal is a model with a single mean density (on the log scale) that is
#  Seasonally dynamic through time, but with a fixed intercept. Sample 
#  variation is Poisson
#

library(pomp)
library(doParallel)


# model components 

  # process model

    process.sim <- "
      Z = b0 + b1 * cos(M_2PI * foy) + b2 * sin(M_2PI * foy);
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
      "

    # from

      fromES <- "
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
               paramnames = c("b0", "b1", "b2", 
                              "Z.0"),
               params = c(b0 = 5, b1 = -0.5, b2 = -0.2,
                          Z.0 = log(Y[1])),
               statenames = c("Z"))


  pf <- pfilter(mod, Np = 1000)
  logLik(pf)



  Z0s <- log(Y[1])
  b0s <- runif(100, 2, 7)
  b1s <- runif(100, -2, 2)
  b2s <- runif(100, -2, 2)
  ep <- data.frame(Z0s, b0s, b1s, b2s)
  colnames(ep) <- c("Z0", "b0", "b1", "b2")


  ncores <- 6
  doParallel::registerDoParallel(ncores)

  lls <- foreach::foreach(i = 1:nrow(ep), .combine = c, 
                           .packages = c("pomp", "magrittr")) %dopar% {
    pfilter(mod, params = c(b0 = ep$b0[i], b1 = ep$b1[i], b2 = ep$b2[i], 
                             Z.0 = ep$Z0[i]), 
             Np = 1000)$loglik

  }

stopImplicitCluster()
par(mfrow = c(4,1))
plot(ep$Z0, lls)
plot(ep$b0, lls)
plot(ep$b1, lls)
plot(ep$b2, lls)






  Z0s <- log(Y[1])
  b0s <- runif(50, 2, 7)
  b1s <- runif(50, -2, 2)
  b2s <- runif(50, -2, 2)
  guesses <- data.frame(Z0s, b0s, b1s, b2s)
  colnames(guesses) <- c("Z0", "b0", "b1", "b2")

  ncores <- 6
  doParallel::registerDoParallel(ncores)

  runs <- foreach::foreach(i = 1:nrow(guesses),  
                           .packages = c("pomp", "magrittr")) %dopar% {
                 mif2(mod, 
                      start = c("Z.0" = guesses[i, 1],
                                "b0" = guesses[i, 2],
                                "b1" = guesses[i, 3],
                                "b2" = guesses[i, 4]),
                      Nmif = 100, Np = 1000, 
                      rw.sd = rw.sd(Z.0 = 0.0001,
                                    b0 = 0.2,
                                    b1 = 0.2,
                                    b2 = 0.2),
                      cooling.fraction.50 = 0.95, 
                      cooling.type = "hyperbolic",
                      transform = TRUE) -> m1
                 m1
  }

stopImplicitCluster()




par(mfrow = c(3, 1))
plot(1, 1, type = 'n', ylim = (c(2, 7)), xlim = c(0, 100), 
     ylab = "b0", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"b0"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}
plot(1, 1, type = 'n', ylim = (c(-3, 3)), xlim = c(0, 100), 
     ylab = "b1", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"b1"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}
plot(1, 1, type = 'n', ylim = (c(-3, 3)), xlim = c(0, 100), 
     ylab = "b2", xlab = "iteration")
for(i in 1:nrow(guesses)){
points((conv.rec(runs[[i]])[,"b2"]), type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}






%>% 
                 continue(Nmif = 200, cooling.fraction = 0.8) %>%
                 continue(Nmif = 200, cooling.fraction = 0.6) %>%
                 continue(Nmif = 200, cooling.fraction = 0.4) %>%
                 continue(Nmif = 200, cooling.fraction= 0.2)  -> m1
                 m1
  }

stopImplicitCluster()




windows(12,12)


