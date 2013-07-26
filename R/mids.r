#'Multiply imputed data set (\code{mids})
#'
#'The \code{mids} object contains a multiply imputed data set. The \code{mids} object is
#'generated by the \code{mice()} and \code{mice.mids()} functions. The \code{mids}
#'class of objects has methods for the following generic functions:
#'\code{print}, \code{summary}, \code{plot}.
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{call}:}{The call that created the object.}
#'    \item{\code{data}:}{A copy of the incomplete data set.}
#'    \item{\code{m}:}{The number of imputations.}
#'    \item{\code{nmis}:}{An array containing the number of missing observations per column.}
#'    \item{\code{imp}:}{A list of \code{ncol(data)} components with the generated multiple
#'imputations. Each part of the list is a \code{nmis[j]} by \code{m} matrix of
#'imputed values for variable \code{j}.}
#'    \item{\code{method}:}{A vector of strings of \code{length(ncol(data))} specifying the
#'elementary imputation method per column.}
#'    \item{\code{predictorMatrix}:}{A square matrix of size \code{ncol(data)}
#'containing integers specifying the predictor set.}
#'    \item{\code{visitSequence}:}{The sequence in which columns are visited.}
#'    \item{\code{form}:}{A vector of strings with length \code{ncol(data)}, specifying
#'formulae. Each string is parsed and executed within the \code{sampler()}
#'function to create terms for the predictor.  An empty string \code{''} means do nothing.
#'The main value
#'lies in the easy specification of interaction terms.}
#'     \item{\code{control}:}{A list with length \code{ncol{data}) with elements \code{NULL} or a
#'list of control parameters for imputation of the corresponnding variable.}
#'    \item{\code{post}:}{A vector of strings of length \code{ncol(data)} with
#'commands for post-processing}
#'    \item{\code{seed}:}{The seed value of the solution.}
#'    \item{\code{iteration}:}{Last Gibbs sampling iteration number.}
#'    \item{\code{lastSeedValue}:}{The most recent seed value.}
#'    \item{\code{loggedEvents}:}{A data.frame with six columns containing a record of
#'automatic removal actions. It is \code{NULL} is no action was made.  At
#'initialization the program does the following three actions:
#'1. A variable that contains missing values, that is not imputed and that is used as a
#'predictor is removed, 2. a constant variable is removed, and 3. a collinear
#'variable is removed. During iteration, the program does the following
#'actions: 1. one or more variables that are linearly dependent are removed
#'(for categorical data, a 'variable' corresponds to a dummy variable), and 2.
#'proportional odds regression imputation that does not converge and is
#'replaced by \code{polyreg}. Column \code{it} is the iteration number at which
#'the record was added, \code{im} is the imputation number, \code{co} is the
#'column number in the data, \code{dep} is the name of the name of the
#'dependent variable, \code{meth} is the imputation method used, and \code{out}
#'is a (possibly long) character vector with the names of the altered or
#'removed predictors.}
#'    \item{\code{extra}:}{A list of lists of extra state or history.  The outer list is indexed
#' by imputed dataset; the inner list by variable.  The contents of the inner list have types
#' given by the specific imputation method of that variable.  Note the indexing scheme differs from that
#' for imp, which can rely on homogenous types across imputation methods.}
#'    \item{\code{chainMean}:}{A list of \code{m} components. Each component is a
#'\code{length(visitSequence)} by \code{maxit} matrix containing the mean of
#'the generated multiple imputations. The array can be used for monitoring
#'convergence.  Note that observed data are not present in this mean.}
#'    \item{\code{chainVar}:}{A list with similar structure of \code{chainMean},
#'containing the covariances of the imputed values.}
#'    \item{\code{prepared}:}{A list containing the input options and data after various cleaning and data checking.
#' That processing may eliminate some degenerate variables from consideration and change the coding of factors.  Normally,
#' this list is only useful for error checking and is only returned if \code{diagnostics} were requested.
#' List members include most of the arguments to \code{\link{mice}}:
#' \code{data}, \code{visitSequence}, \code{method}, \code{defaultMethod}, \code{predictorMatrix}, \code{form},
#' \code{control}, and \code{post}.  Note that \code{form} will be a list of formula objects, present for all
#' variables to be imputed, rather than the character vector that it is on input.
#' The list also includes
#' \code{nmis} (number missing in each column of data), \code{nvar} (number of columns of data),
#' \code{varnames} (variable names for the data).}
#'\item{\code{loggedEvents}:}{As above.}
#'}
#'
#' @note Many of the functions of the \code{mice} package do not use the S4 class definitions,
#' and instead rely on the S3 list equivalent \code{oldClass(obj) <- "mids"}.
#'
#' The duplication of information between the \code{mids} object and its \code{prepared} member is undesirable.
#' While there may be some differences between the two, they are much less than with the padded (with dummies)
#' data that was the historical basis of this design.
#'
#'@name mids-class
#'@rdname mids-class
#'@aliases mids-class mids
#'@author Stef van Buuren, Karin Groothuis-Oudshoorn, 2000
#'@seealso \code{\link{mice}}, \code{\link[=mira-class]{mira}}, \code{\link[=mipo-class]{mipo}}
#'@references van Buuren S and Groothuis-Oudshoorn K (2011). \code{mice}:
#'Multivariate Imputation by Chained Equations in \code{R}. \emph{Journal of
#'Statistical Software}, \bold{45}(3), 1-67.
#'\url{http://www.jstatsoft.org/v45/i03/}
#'@keywords classes
#'@export
setClass("mids",
         representation(
             call      = "call",
             data      = "data.frame" ,
             m         = "numeric",
             nmis      = "integer",
             imp       = "list",
             method    = "character",
             predictorMatrix = "matrix",
             visitSequence = "numeric",
             form      = "character",
             control   = "list",
             post      = "character",
             seed      = "numeric",
             iteration = "integer",
             lastSeedValue = "integer",
             loggedEvents = "data.frame",
             extra     = "list",
             chainMean = "array",
             chainVar  = "array",
             prepared  = "list"),
         contains  = "list"
)
