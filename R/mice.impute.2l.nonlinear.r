require("plyr")

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

mice.impute.2lmixed.logit <- function(y, ry, x, type, intercept=TRUE, ...)
{
    ## mixed level 1 and 2 predictors of outcomes
    ## variables with 2 at end of name are level 2, one obs /cluster

    ## We make two additions to the input data.
    ## We assume there is an unobserved variable z with
    ## Pr(y=1) = 1/(1+exp(-z))
    ## Second, for each group there is a continuous group effect theta2
    ## with mean 0 and sd tau
    ## z=Xb+theta+e
    ## (theta is theta2 duplicated for each observation at level 1)
    ## e is normal error with sd sigma
    ## It seems likely that sigma and tau may be correlated, slowing convergence.
  ## Initialize
  n.iter <- 5
  #n.iter <- 1000
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
  sigma <- NA_real_
  z.prior.mean = rep(NA_real_, nrow(X))


# for each iteration record sigma2, beta.post.mean, beta, mu2, y2
  nvar <- ncol(v)
  ntrace <- 2+nvar+n.class+nrow(X)+nmiss+nvar+nrow(X)
  MCTRACE <<- matrix(NA_real_, nrow=n.iter+1, ncol= ntrace)
  # order is all the posterior values and then some related stats
  MCTRACE[1,] <<- c(tau, sigma, beta, theta2, z, y[nry], beta.post.mean, z.prior.mean)

  #optimization: precompute constant for inner loop posterior()
  grid.lo <- -3
  grid.hi <- 3
  grid.size <- 250
  grid.raw <- seq(grid.lo, grid.hi, length=grid.size)
  # The true grid is mu+sigma*grid.raw and we need the probabilities
  # at those points.  But the probability for N(mu, sd) at mu+sd*x
  # is [probability for N(0, 1) at x]/sd.  The division
  # is just an additive constant for log(prob) and can be ignored.
  # So we only need to compute the probabilities once. dnorm is relatively expensive.
  grid.lnprob <- log(dnorm(grid.raw))

  for (iter in 1:n.iter){
      # X already has an intercept in it

      # Gelman et al 2004 pp. 299-301 for hierarchical normal model
      # draw posterior values for tau given all other values
      # theta2 ~ Norm(0, tau)
      tau <- sqrt(sum(theta2^2)/rchisq(1, n.class-1))

      # and for sigma, also known mean 0
      resid <- z- X %*% beta - theta
      sigma <- sqrt(sum(resid^2)/rchisq(1, nrow(X)))
      # the use of df here and df-1 above is deliberate
      # reflecting different prior distns needed for proper posterior

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

      ## draw posterior beta and redraw sigma
      w <- z-theta
      beta.post.mean <- t(w %*% X %*% v)
      residuals <- z - X %*% beta.post.mean - theta
      sigma <- sqrt(sum((residuals)^2)/rchisq(1, nrow(X) - ncol(X)))  # SvB 01/02/2011
      beta <- beta.post.mean + (t(chol((v + t(v))/2)) %*% rnorm(ncol(X))) * sigma

      # z.prior.mean is center of prior dist of z given the coefficient draws
      z.prior.mean <- X %*% beta + theta

      # Prior z ~ N(z.prior.mean, sigma)
      # commence computation of posterior (z|y and all other values)
      # done numerically

      posterior <- function(x){
          mu <- x[1]  #["z.prior.mean"] except cbind does not label column
          sigma <- x["sigma"]
          y <- x["y"]
          grid <- mu + sigma*grid.raw
          p = 1/(1+exp(-grid))
          if ( y == 0)
              p = 1- p
          post <- log(p) + grid.lnprob
          post <- post-max(post)
          sample(grid, 1, prob=exp(post))
      }
      z <- apply(cbind(z.prior.mean, sigma, y), 1, posterior)
      # and impute the missing observed values
      y[nry] <- rbinom(nmiss, 1, 1/(1+exp(-z[nry])))
      MCTRACE[iter+1,] <<- c(tau, sigma, beta, theta2, z, y[nry], beta.post.mean, z.prior.mean)
  }
  return(y[!ry])
}
