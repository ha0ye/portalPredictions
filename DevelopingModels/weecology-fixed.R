#
# Weecology-Fixed
#
# Fixed is a model with a single mean density (on the log scale) that does
#  not evolve over time. Sample variation is Poisson
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

      Y = rpois(exp(Z));
    "

    count.dmeas <- "

      lik = dpois(Y, exp(Z), give_log);
    "

  mod <- pomp(data = data.frame(time = nm, Y = Y), 
               times = "time", t0 = 216,
               rprocess = discrete.time.sim(step.fun = Csnippet(process.sim),
                                            delta.t = 1), 
               rmeasure = Csnippet(count.rmeas),
               dmeasure = Csnippet(count.dmeas),
               paramnames = c("Z.0"),
               params = c(Z.0 = log(mean(na.omit(Y)))),
               statenames = c("Z"))


  pf <- pfilter(mod, Np = 1000)
  logLik(pf)

  Z0s <- seq(2, 7, 0.1)

  ncores <- 6
  doParallel::registerDoParallel(ncores)

  lls <- foreach::foreach(i = 1:length(Z0s), .combine = rbind, 
                           .packages = c("pomp", "magrittr")) %dopar% {

     pfilter(mod, params = c(Z.0 = Z0s[i]), Np = 1000)$loglik
  }

lls
plot(Z0s, lls[,1])


  guesses <- data.frame(Z.0 = seq(2, 7, 0.5))

 


  ncores <- 6
  doParallel::registerDoParallel(ncores)

  runs <- foreach::foreach(i = 1:nrow(guesses), .inorder = FALSE, 
                           .packages = c("pomp", "magrittr")) %dopar% {

                 mif2(mod, start = c("Z.0" = guesses[i, 1]),
                      Nmif = 50, Np = 1000, rw.sd = rw.sd(Z.0 = 0.5),
                      cooling.fraction.50 = 0.95, 
                      cooling.type = "hyperbolic") %>% 
                 continue(Nmif = 50, cooling.fraction = 0.8) %>%
                 continue(Nmif = 50, cooling.fraction = 0.6) %>%
                 continue(Nmif = 50, cooling.fraction= 0.2) -> m1
                 m1
  }


windows(12,12)
plot(1, 1, type = 'n', ylim = c(3.8, 5), xlim = c(0, 1000))
for(i in 1:nrow(guesses)){
points(conv.rec(runs[[i]])[,"Z.0"], type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}


for(i in 1:nrow(guesses)){
                 mif2(mod, start = c("Z.0" = guesses[i, 1]),
                      Nmif = 100, Np = 1000, rw.sd = rw.sd(Z.0 = 0.01),
                      cooling.fraction.50 = 0.6, 
                      cooling.type = "hyperbolic") %>% 
                 mif2() -> mf
                 runs[[i]] <- mf
}


windows(12,12)
plot(1, 1, type = 'n', ylim = c(2,6), xlim = c(0, 100))
for(i in 1:nrow(guesses)){
points(conv.rec(runs[[i]])[,"Z.0"], type= 'l', 
       lwd = 2, col = rgb(0.6, 0.6, 0.6, 0.6))
}




  runs <- foreach::foreach(i = 1:nrow(guesses), 
                           .combine = rbind, 
                           .packages = c("pomp", "magrittr")) %dopar% {

                 mif2(mod, start = c("Z.0" = guesses[i, 1]),
                      Nmif = 50, Np = 1000, rw.sd = rw.sd(Z.0 = 0.01),
                      cooling.fraction.50 = 0.8, 
                      cooling.type = "geometric") %>% 
                 mif2() -> mf
                 ll <- logmeanexp(replicate(5, logLik(pfilter(mf))), 
                                  se = TRUE)
                 data.frame(loglik = ll[1], loglik.se = ll[2],
                            as.list(coef(mf)))
  }
    





