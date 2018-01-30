library(pomp)
library(magrittr)
library(foreach)
library(doParallel)

pompExample(gompertz)

plot(gompertz)

theta <- coef(gompertz)
theta.true <- theta

registerDoParallel()

estpars <- c("r", "sigma", "tau")

foreach(i=1:10,
        .inorder=FALSE,
        .packages = c("magrittr", "pomp"),
        .options.multicore=list(set.seed=TRUE)
        ) %dopar%
    {
        theta.guess <- theta.true
        theta.guess[estpars] <- rlnorm(
            n=length(estpars),
            meanlog=log(theta.guess[estpars]),
            sdlog=1
        )
        mif2(
            gompertz,
            Nmif=50,
            start=theta.guess,
            transform=TRUE,
            rw.sd=rw.sd(r=0.02,sigma=0.02,tau=0.05),
            cooling.fraction.50=0.95,
            Np=2000
        ) %>%
            continue(Nmif=50,cooling.fraction=0.8) %>%
            continue(Nmif=50,cooling.fraction=0.6) %>%
            continue(Nmif=50,cooling.fraction=0.2) -> m1
        ll <- replicate(n=10,logLik(pfilter(m1,Np=10000)))
        list(mif=m1,ll=logmeanexp(ll,se=TRUE))
    } -> mf