\name{rbind.mids}
\alias{rbind.mids}
\title{Rowwise combination of \code{mids} objects from runs on partial data.}
\usage{
  rbind.mids(x, y, ...)
}
\arguments{
  \item{x}{A \code{mids} object.}

  \item{y}{A \code{mids} object or a \code{data.frame},
  \code{matrix}, \code{factor} or \code{vector}.}

  \item{\dots}{Additional \code{data.frame}, \code{matrix},
  \code{vector} or \code{factor}. These can be given as
  named arguments.}
}
\value{
  An S3 object of class \code{mids} Component \code{call}
  is a vector, with first argument the \code{mice()}
  statement that created \code{x} and second argument the
  call to \code{rbind.mids()}. Component \code{data} is the
  rowwise combination of the (incomplete) data in \code{x}
  and \code{y}. Component \code{nmis} is an array
  containing the number of missing observations per column,
  defined as \code{x$nmis} + \code{y$nmis}. Component
  \code{imp} is a list of \code{nvar} components with the
  generated multiple imputations.  Each part of the list is
  a \code{nmis[j]} by \code{m} matrix of imputed values for
  variable \code{j}. If \code{y} is a \code{mids} object
  then \code{imp[[j]]} equals \code{rbind(x$imp[[j]],
  y$imp[[j]])}; otherwise the original data of \code{y}
  will be copied into this list, including the missing
  values of \code{y} then \code{y} is not imputed.
  Component \code{chainMean} is set to \code{NA} if
  \code{y} is \code{mids}, otherwise \code{x$chainMean}
  (which is probably not reliable). Component
  \code{chainVar} is set to \code{NA} if \code{y} is
  \code{mids}, otherwise \code{x$chainMean} (which is
  probably not reliable). Component \code{prepared} is set
  to \code{x$prepared} unless \code{y} is \code{mids}, in
  which case the \code{prepared$data} elements are
  combined. Component \code{loggedEvents} is set to
  \code{x$loggedEvents}. Other elements are simply copied
  from \code{x}.
}
\description{
  Append \code{mids} objects by rows, assuming they
  represent 2 runs on different subsets of the data.
}
\details{
  If both \code{x} and \code{y} are \code{mids} objects the
  function requires that they represent runs with the same
  specification (arguments to \code{mice}), except that
  they are run on different subsets (by rows) of the data.
  The resulting \code{mids} object combines the data and
  imputations from both runs, and either combines other
  information or copies it from the first argument.

  If \code{y) is not a \code{mids} object, \code{x$data},
  \code{y} and the optional arguments \code{...} are
  stacked with \code{\link{rbind}} and treated as the new
  data.  Imputations will be missing for all rows except
  those corresponding to \code{x}.
}
\note{
  The results of imputations on subsets of the data, even
  when combined, will not be the same as the result of
  imputation on the entire data.
}
\author{
  Karin Groothuis-Oudshoorn, Stef van Buuren, 2009
}
\references{
  van Buuren S and Groothuis-Oudshoorn K (2011).
  \code{mice}: Multivariate Imputation by Chained Equations
  in \code{R}. \emph{Journal of Statistical Software},
  \bold{45}(3), 1-67.
  \url{http://www.jstatsoft.org/v45/i03/}
}
\seealso{
  \code{\link{cbind.mids}}, \code{\link{ibind}},
  \code{\link[=mids-class]{mids}}
}
\keyword{manip}

