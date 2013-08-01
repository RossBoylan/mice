\name{mice.impute.2l.logit}
\alias{2l.logit}
\alias{2l.probit}
\alias{mice.impute.2l.logit}
\alias{mice.impute.2l.probit}
\title{Imputation of binary outcome by two-level model. Not necessarily logistic!}
\usage{
  mice.impute.2l.logit(y, ry, x, type, intercept = TRUE,
    ...)
}
\arguments{
  \item{y}{Incomplete data vector of length \code{n} of 0's
  and 1's.}

  \item{ry}{Vector of missing data pattern
  (\code{FALSE}=missing, \code{TRUE}=observed)}

  \item{x}{Matrix (\code{n} x \code{p}) of complete
  covariates.}

  \item{type}{Vector of length \code{ncol(x)} identifying
  random and class variables.  Random variables are
  identified by a '2'. The class variable (only one is
  allowed) is coded as '-2'. Random variables also include
  the fixed effect.}

  \item{intercept}{Logical determining whether the
  intercept is automatically added.}

  \item{...}{Other named arguments.}
}
\value{
  A vector of length \code{nmis} with imputations.
}
\description{
  Imputes a binary outcome, coded 0/1, using a two-level
  model
}
\details{
  If all the covariates are at level 2 (group) this does
  logistic regression. In the more typical case in which
  some of the variables are level 1 this currently uses a
  probit model with Hybrid Monte-Carlo. Code for the Albert
  and Chib estimator using latent Z is also provided,
  though it is not currently hooked up.
}
\author{
  Ross D. Boylan
}

