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
  #n.iter <- 100
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

# mixed level 1 and 2 predictors of outcomes
mice.impute.2lmixed.logit <- function(y, ry, x, type, intercept=TRUE, ...)
{
    ## We make two additions to the input data.
    ## We assume there is an unobserved variable z with
    ## Pr(y=1) = 1/(1+exp(-z))
    ## Second, for each group there is a continuous group effect z2
    ## with mean 0 and variance tauz
    ## It seems likely that tau and tauz may be correlated, slowing convergence.
    ## z=Xb+z2
    ## Define zres= z-Xb
  ## Initialize
  #n.iter <- 100
  n.iter <- 1000
  nry <- !ry
  n.class <- length(unique(x[, type==(-2)]))
  gf.full <- factor(x[,type==(-2)], labels=1:n.class)
  ids <- contract(data.frame(gf.full), gf.full)
  gf <- gf.full[ry]

  X <- as.matrix(x[,type>0])

  # level 1 latent variables
  # establish initial guesses with Laplace's law of succession
  p <- (y+1)/3
  # set missing values to the sample average
  p[is.na(p)] <- mean(y, na.rm=TRUE)
  # invert to get the latent variable
  z <- qlogis(p)
  zres <- z-mean(z)  # inital guess assumes all coefficient except intercept are 0

  # at this point we have imputed values for the latent z
  # but not for the y

  # compute some constants for the loop
  xtx <- t(X) %*% X
  ridge <- 0.00001
  pen <- ridge * diag(xtx)
  if (length(pen)==1) pen <- matrix(pen)
  v <- solve(xtx + diag(pen))

  # continue setting initial values
  coef <- t(z %*% X %*% v)
  # for each iteration record sigma2, coef, beta, mu2, y2
  nvar <- ncol(v)
  ntrace <- 1+nvar+nvar+n.class+n.class
  MCTRACE <<- matrix(NA_real_, nrow=n.iter+1, ncol= ntrace)
  MCTRACE[1, seq(ntrace-n.class+1, ntrace)] <<- y2

  for (iter in 1:n.iter){
      # X already has an intercept in it

      # WORK IN PROGRESS

      # draw values for z2 given all other values
      theta2 <- c(by(y, gf.full, mean))
      theta2e <- expand(theta2, gf.full,

      # the next section is modeled on .norm.draw except
      # it estimates the model from all the data
      # and then imputes all the data.
      # Also, this works on the level 2 data.

      coef <- t(y2 %*% X %*% v)
      residuals <- y2 - X %*% coef
      sigma2 <- sqrt(sum((residuals)^2)/rchisq(1, nrow(X) - ncol(X)))  # SvB 01/02/2011
      beta <- coef + (t(chol((v + t(v))/2)) %*% rnorm(ncol(X))) * sigma2

      # mu2 is center of prior dist of y2 given the coefficient draws
      mu2 <- X %*% beta

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
