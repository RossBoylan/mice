require("plyr")
require("truncnorm")
require("HybridMC")

## exploratory functions to get posterior likelihoods given a linear predictor
## that is normal
## y is a vector of 0/1 outcomes.
# theta is a vector of points at which to evaluate the posterior
# note the result is not normalized and is a log likelihood
binpost <- function(theta, y, mu, sigma){
    p <- 1/(1+exp(-theta))
    nyes <- sum(y, na.rm=TRUE)
    (nyes*log(p)+(sum(!is.na(y))-nyes)*log(1-p))+log(dnorm(theta, mean=mu, sd=sigma))
}

# return only one row of data per factor level
# uses first available row
# f is a factor with length equal to the number of rows in data
contract <- function(data, f){
    i <- match(levels(f), f)
    if (is.null(dim(data)))
        data.frame(data[i])
    else
        data[i,]
}

# expanded has multiple rows per factor
# contracted has one row per factor
# f[i] is  the factor for expanded[i,]
# replicate the info in contracted into expanded
# Assumes factors for contracted are in same order as for f
# both is an optional argument giving column names to copy.
# They must be present in both datasets
expand <- function(contracted, expanded, f,
                   both=intersect(colnames(contracted), colnames(expanded))){
    i <- match(levels(f), f)
    cf <- f[i]  # factors in contracted
    expanded[,both] <- contracted[match(f, cf), both]
    expanded
}


#-------------------MICE.IMPUTE.2l.logit----------------------------

mice.impute.2l.logit <- function(y, ry, x, type, intercept=TRUE, ...)
{

  ## append intercept
  if (intercept) {
    x <- cbind(1, as.matrix(x))
    type <- c(2, type)
  }

  if (any(type==1))
      return (mice.impute.2lmixed.logit(y, ry, x, type))


  ## Initialize
  #n.iter <- 5
  n.iter <- 1000
  nry <- !ry
  n.class <- length(unique(x[, type==(-2)]))
  gf.full <- factor(x[,type==(-2)], labels=1:n.class)
  gf <- gf.full[ry]

  # 2 is a suffix for 2nd level, i.e., group-level, material
  X2 <- as.matrix(contract(x[,type==2], gf.full))
  id2 <- contract(data.frame(gf.full), gf.full)
  r2 <- rep(FALSE, nrow(X2))  # latent variable, never observed
  p2 <- c(by(y, gf.full, function(x) mean(x, na.rm=TRUE)))
  p2[is.na(p2)] <- mean(y, na.rm=TRUE)
  p2[p2>0.95] <- 0.95
  p2[p2<0.05] <- 0.05
  # y2 is the latent linear predictor in a logistic model
  # note that it is continuous, while y is binary
  y2 <- qlogis(p2)

  # at this point we have imputed values for all the y2
  # but not for the y

  # compute some constants for the loop
      xtx <- t(X2) %*% X2
      ridge <- 0.00001
      pen <- ridge * diag(xtx)
      if (length(pen)==1) pen <- matrix(pen)
      v <- solve(xtx + diag(pen))

  # for each iteration record sigma2, coef, beta, mu2, y2
  nvar <- ncol(v)
  ntrace <- 1+nvar+nvar+n.class+n.class
  MCTRACE <<- matrix(NA_real_, nrow=n.iter+1, ncol= ntrace)
  MCTRACE[1, seq(ntrace-n.class+1, ntrace)] <<- y2

  for (iter in 1:n.iter){
      # X2 already has an intercept in it

      # the next section is modeled on .norm.draw except
      # it estimates the model from all the data
      # and then imputes all the data.
      # Also, this works on the level 2 data.

      coef <- t(y2 %*% X2 %*% v)
      residuals <- y2 - X2 %*% coef
      sigma2 <- sqrt(sum((residuals)^2)/rchisq(1, nrow(X2) - ncol(X2)))  # SvB 01/02/2011
      beta <- coef + (t(chol((v + t(v))/2)) %*% rnorm(ncol(X2))) * sigma2

      # mu2 is center of prior dist of y2 given the coefficient draws
      mu2 <- X2 %*% beta

      # Prior y2 ~ N(mu2, sigma2)  Note sigma2 is 2nd level, not squared
      # We do not consider the uncertainty in mu2 and sigma2
      # commence computation of posterior y2
      # TODO: stop using column indices to refer to values in the argument of posterior
      count2 <- t(simplify2array(by(y, gf.full, function(b) c(sum(b, na.rm=TRUE), sum(!is.na(b))))))
      # draw a posterior y for each group
      # Numerically intensive
      posterior <- function(x) {
          nyes <- x[1]
          nno <- x[2]-nyes
          mu <- x[3]
          grid <- seq(mu-3*sigma2, mu+3*sigma2, length.out=500)
          p <- 1/(1+exp(-grid))
          post <- (nyes*log(p)+nno*log(1-p))+log(dnorm(grid, mean=mu, sd=sigma2))
          post <- post-max(post)
          sample(grid, 1, prob=exp(post))
      }
      # impute all the  posterior ys
      y2 <- apply(cbind(count2, mu2), 1, posterior)
      p2 <- 1/(1+exp(-y2))

      # apply p2 to impute y for each missing value
      # these are indices into p2 from all the missing observations in y
      into2 <- match(gf.full[nry], id2)
      y[nry] <- rbinom(sum(nry), 1, p2[into2])
                                        # for each iteration record sigma2, coef, beta, mu2, y2
      MCTRACE[iter+1,] <<- c(sigma2, coef, beta, mu2, y2)
  }
  return(y[!ry])
}

mice.impute.2lmixed.logit.AlbertChib <- function(y, ry, x, type, intercept=TRUE, ...)
{
    ## mixed level 1 and 2 predictors of outcomes
    ## variables with 2 at end of name are level 2, one obs /cluster

    ## We make two additions to the input data.
    ## We assume there is an unobserved variable z with
    ## y=1 iff z>0.
    ## This specification is from Albert and Chib 1993, though
    ## that model was single level only.  We extend it ...
    ## Second, for each group there is a continuous group effect theta2
    ## with mean 0 and sd tau
    ## z=Xb+theta
    ## (theta is theta2 duplicated for each observation at level 1)

  ## Initialize
  #n.iter <- 5
  n.iter <- 1000
  nry <- !ry
  nmiss <- sum(nry)
  n.class <- length(unique(x[, type==(-2)]))
  gf.full <- factor(x[,type==(-2)], labels=1:n.class)
  ids <- contract(data.frame(gf.full), gf.full)
  gf <- gf.full[ry]

  X <- as.matrix(x[,type>0])

  # compute some constants for the loop
  xtx <- t(X) %*% X
  ridge <- 0.00001
  pen <- ridge * diag(xtx)
  if (length(pen)==1) pen <- matrix(pen)
  v <- solve(xtx + diag(pen))

  # level 1 latent variables
  # establish initial guesses with Laplace's law of succession
  p <- (y+1)/3
  # set missing values to the sample average
  # I think that y have already been imputed before entering this function
  p[is.na(p)] <- mean(y, na.rm=TRUE)
  # invert to get the latent variable
  z <- qlogis(p)

  # at this point we have imputed values for the latent z
  # but not for the y.  I think it's unnecessary, but just in case...
  y[nry] <- sample(y[ry], nmiss, replace=TRUE)

  # continue setting initial values

  beta.post.mean <- t(z %*% X %*% v)
  beta <- beta.post.mean
  resid <- z-X%*%beta
  # level 2 latent variables initial values
  theta2 <- ddply(data.frame(id=gf.full, resid=resid), .(id), summarize, mean=mean(resid))
  iExpand <- match(gf.full, theta2$id)
  theta2 <- theta2[,"mean"]
  theta <- theta2[iExpand]

  # No initial values needed
  tau <- NA_real_
  z.prior.mean = rep(NA_real_, nrow(X))

  # level 1 variance is a constant in this model
  sigma <- 1.0

  nvar <- ncol(v)
  ntrace <- 1+nvar+n.class+nrow(X)+nmiss+nvar+nrow(X)
  MCTRACE <<- matrix(NA_real_, nrow=n.iter+1, ncol= ntrace)
  # order is all the posterior values and then some related stats
  # final value is acceptance rate for candidate z
  MCTRACE[1,] <<- c(tau, beta, theta2, z, y[nry], beta.post.mean, z.prior.mean)

  # pull calculations out of loop
  # if y is missing we draw a normal, otherwise a truncated normal
  # consistent with the observed y
  trunclo <- ifelse(ry & (y==1), 0, -Inf)
  trunchi <- ifelse(ry & (y==0), 0, Inf)

  for (iter in 1:n.iter){
      # X already has an intercept in it

      # Gelman et al 2004 pp. 299-301 for hierarchical normal model
      # draw posterior values for tau given all other values
      # theta2 ~ Norm(0, tau)
      tau <- sqrt(sum(theta2^2)/rchisq(1, n.class-1))

      # draw posterior theta
      thetapost <- ddply(data.frame(id=gf.full, r= z - X %*%beta.post.mean), .(id), function(df) {
          x <- df$r
          n <- length(x)
          varpost <- 1/(1/tau^2 + n/sigma^2)
          mupost <- sum(x)/sigma^2*varpost
          data.frame(mu=mupost, var=varpost)
      })
      theta2 <- rnorm(n.class, mean=thetapost$mu,
                      sd=sqrt(thetapost$var))
      # replicate theta2 onto level 1
      theta <- theta2[iExpand]

      ## draw posterior beta
      w <- z-theta
      beta.post.mean <- t(w %*% X %*% v)
      residuals <- z - X %*% beta.post.mean - theta
      beta <- beta.post.mean + (t(chol((v + t(v))/2)) %*% rnorm(ncol(X))) * sigma

      # z.prior.mean is mean parameter for center of prior dist of z given the coefficient draws
      z.prior.mean <- X %*% beta + theta
      z <- rtruncnorm(nrow(X), a=trunclo, b=trunchi, mean=z.prior.mean)

      # and impute the missing observed values
      y[nry] <- z[nry]>0
      MCTRACE[iter+1,] <<- c(tau, beta, theta2, z, y[nry], beta.post.mean, z.prior.mean)
  }
  return(y[!ry])
}

# functions for use by hybrid monte carlo
# These are the log density of a probit model (without any Z terms)
# and its derivatives
# parameter ordering for q is beta, tau, theta2
logDens <- function(q, nvar, n.class, gf.full, X, y, iExpand) {
    beta <- q[1:nvar]
    theta2 <- q[(nvar+2):length(q)]
    tau <- q[nvar+1]
    eta <- X %*% beta + theta2[iExpand]
    lprob <- ifelse(y, pnorm(eta, lower.tail=TRUE, log.p=TRUE),
                    pnorm(eta, lower.tail=FALSE, log.p=TRUE))
    f <- sum(lprob)+sum(dnorm(theta2, sd=tau, log=TRUE))-2*log(tau)
    f
}

dLogDens <- function(q, nvar, n.class, gf.full, X, y, iExpand) {
      if(any(is.na(q))) recover()
      beta <- q[1:nvar]
      theta2 <- q[(nvar+2):length(q)]
      tau <- q[nvar+1]
      theta <- theta2[iExpand]
                                        # compute derivatives using formula
      eta <- X %*% beta + theta
      cum <- pnorm(eta)
      dens <- dnorm(eta)
      # we want term1 = (y-cum)*dens/[cum*(1-cum)]
      # this can get ugly at extreme eta
      # if y=1, term1 = dens/cum (ugly for eta<<0)
      # if y=0, term1 = -dens/(1-cum) (ugly for eta>>0)
      # in either ugly case, Hopital + maxima says it behaves roughly like abs(eta)
      iWild <- ifelse(y==1, eta < -4, eta > 4)
      term1 <- abs(eta)  # easiest to set dimensions with the extreme case
      iRegular <- (! iWild ) & (y==1)
      # subset at start to minimize computations
      term1[iRegular] <- dens[iRegular]/cum[iRegular]

      iRegular <- (! iWild) & (y==0)
      term1[iRegular] <- dens[iRegular]/(cum[iRegular]-1)

      term1mat <- matrix(term1, nrow=nrow(X), ncol=nvar)
      dldbeta <- apply(term1mat*X, 2, sum)
      # tapply and contract, which orders theta2, both
      # use the same order, and so the addition below should work.
      dldtheta <- tapply(term1, gf.full, sum)-theta2/tau^2
      dldtau <- (-2-n.class+sum(theta2^2)/tau^2)/tau
      r <- c(dldbeta, dldtau, dldtheta)
    if(any(is.na(r))) recover()
      MCTRACE2 <<- c(MCTRACE2, list(c(r, q)))
      r
      #r <- pmax(-500, r)
      #pmin(500, r)
  }

# calculate numeric derivative of scalar value f at x
# by varying the i'th component of x by steps
slope <- function(f, x, i, steps) {
    delta <- rep(0, length(x))
    r0 <- f(x)
    r1 <- sapply(steps, function(s) {delta[i] <- delta[i]+s; f(x+delta)})
    bump <- r1-r0
    data.frame(step=c(0, steps), f=c(r0, r1), deriv=c(NA, bump/steps))
}

# Use hybrid monte carlo
# Initially I'll just see if I can compute the derivatives correctly
mice.impute.2lmixed.logit <- function(y, ry, x, type, intercept=TRUE, ...)
{
    ## mixed level 1 and 2 predictors of outcomes
    ## variables with 2 at end of name are level 2, one obs /cluster


  ## Initialize
  MCTRACE2 <<- list()  # will hold calls to logDens
  #n.iter <- 5
  n.iter <- 1000
  nry <- !ry
  nmiss <- sum(nry)
  n.class <- length(unique(x[, type==(-2)]))
  gf.full <- factor(x[,type==(-2)], labels=1:n.class)
  ids <- contract(data.frame(gf.full), gf.full)
  iExpand <- match(gf.full, ids)
  gf <- gf.full[ry]

  X <- as.matrix(x[,type>0])
  nvar <- ncol(X)

  # pick plausible starting values by doing a linear regression
  # on a logistic transform of observation
  xtx <- t(X) %*% X
  ridge <- 0.00001
  pen <- ridge * diag(xtx)
  if (length(pen)==1) pen <- matrix(pen)
  v <- solve(xtx + diag(pen))

  # level 1 latent variables
  # establish initial guesses with Laplace's law of succession
  p <- (y+1)/3
  # set missing values to the sample average
  # I think that y have already been imputed before entering this function
  p[is.na(p)] <- mean(y, na.rm=TRUE)
  # invert to get the latent variable
  z <- qlogis(p)

  # at this point we have imputed values for the latent z
  # but not for the y.  I think it's unnecessary, but just in case...
  y[nry] <- sample(y[ry], nmiss, replace=TRUE)

  # continue setting initial values

  beta <- t(z %*% X %*% v)
  resid <- z-X%*%beta
  # level 2 latent variables initial values
  theta2 <- ddply(data.frame(id=gf.full, resid=resid), .(id), summarize, mean=mean(resid))
  theta2 <- theta2[,"mean"]
  tau <- sqrt(sum((theta2-mean(theta2))^2)/(n.class-nvar))

  ## # check analytic derivatives
  ## danalytic <- dLogDens(c(beta, tau, theta2), nvar, n.class, gf.full, X, y, iExpand)
  ## myenv <- new.env()
  ## assign("x", c(beta, tau, theta2), envir=myenv)
  ## assign("nvar", nvar, envir=myenv)
  ## assign("n.class", n.class, envir=myenv)
  ## assign("gf.full", gf.full, envir=myenv)
  ## assign("X", X, envir=myenv)
  ## assign("y", y, envir=myenv)
  ## assign("iExpand", iExpand, envir=myenv)
  ## assign("logDens", logDens, envir=myenv)
  ## r <- numericDeriv(quote(logDens(x, nvar, n.class, gf.full, X, y, iExpand)), "x", myenv)
  ## dnumeric <- c(attr(r, "gradient"))
  ## delta <- danalytic-dnumeric
  ## tempfun <- function(x) logDens(x, nvar, n.class, gf.full, X, y, iExpand)
  ## myderiv <- slope(tempfun, c(beta, tau, theta2), 1, 10^seq(-5, 2))
  ## browser()

  # scale the problem
  # didn't seem to work, ie. chain didn't move.  got lots of warnings that
  #49: In dnorm(theta2, sd = tau, log = TRUE) : NaNs produced
  #50: In log(tau) : NaNs produced

  q <- c(beta, tau, theta2)
  danalytic <- dLogDens(q, nvar, n.class, gf.full, X, y, iExpand)
  myenv <- new.env()
  assign("x", q, envir=myenv)
  assign("nvar", nvar, envir=myenv)
  assign("n.class", n.class, envir=myenv)
  assign("gf.full", gf.full, envir=myenv)
  assign("X", X, envir=myenv)
  assign("y", y, envir=myenv)
  assign("iExpand", iExpand, envir=myenv)
  assign("dLogDens", dLogDens, envir=myenv)
  # the next step is slow since it computes all cross derivatives
  r <- numericDeriv(quote(dLogDens(x, nvar, n.class, gf.full, X, y, iExpand)), "x", myenv)
  d2 <- diag(matrix(c(attr(r, "gradient")), nrow=length(q)))
  iZero <- d2 == 0
  d2a <- ifelse(iZero, danalytic/min(abs(d2[!iZero]))/2, danalytic/d2)
  weights <- abs(1/d2a)
  #weights <- d2a^2
  #weights <- 1.0
  
  #epsilon <- c(0.01, 0.04)/2
  epsilon <- 0.01
  LFsteps <- 30
  r <- HybridMC::hybridMC(y.start=c(beta, tau, theta2), n.samp=100,
                          logDens=logDens, dLogDens=dLogDens, epsilon=epsilon,
                          LFsteps=LFsteps, compWeights=weights, MPwidth=1,
                          MPweights=1,
                          nvar=nvar, n.class=n.class, gf.full=gf.full,
                          X=X, y=y,
                          iExpand=iExpand)
  # note this is already type mcmc
  MCTRACE <<- r

  # impute missing y from final value of parameters
  q <- r[end(r),]
  beta <- q[1:nvar]
  theta2 <- p[(nvar+1):length(q)]
  tau <- p[nvar+1]
  theta <- theta2[iExpand]

  eta <- X[!ry,] %*% beta + theta[!ry]
  ymiss <- rbinom(nmiss, 1, pnorm(eta))
  MCTRACE2 <<- do.call(rbind, MCTRACE2)
  return(ymiss)
}
