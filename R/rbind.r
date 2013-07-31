#'Rowwise combination of \code{mids} objects from runs on partial data.
#'
#'Append \code{mids} objects by rows, assuming they represent 2 runs on different
#' subsets of the data.
#'
#'If both \code{x} and \code{y} are \code{mids} objects the function requires that they
#' represent runs  with the same specification (arguments to \code{mice}), except that
#' they are run on different subsets (by rows) of the data.  The resulting \code{mids} object combines
#' the data and imputations from both runs, and either combines other information
#' or copies it from the first argument.
#'
#' If \code{y} is not a \code{mids} object, \code{x$data}, \code{y} and the optional arguments
#' \code{...} are stacked with \code{\link{rbind}} and treated as the new data.  Imputations
#' will be missing for all rows except those corresponding to \code{x}.
#'@param x A \code{mids} object.
#'@param y A \code{mids} object or a \code{data.frame}, \code{matrix}, \code{factor}
#'or \code{vector}.
#'@param \dots Additional \code{data.frame}, \code{matrix}, \code{vector} or \code{factor}.
#'These can be given as named arguments.
#'@return An S3 object of class \code{mids}
#'Component \code{call} is a vector, with first argument the \code{mice()} statement
#'that created \code{x} and second argument the call to \code{rbind.mids()}. Component
#'\code{data} is the rowwise combination of the (incomplete) data in \code{x}
#'and \code{y}.
#'Component \code{nmis} is an array containing the number of missing observations per
#'column, defined as \code{x$nmis} + \code{y$nmis}.
#'Component \code{imp} is a list of \code{nvar} components with the generated multiple
#'imputations.  Each part of the list is a \code{nmis[j]} by \code{m} matrix of
#'imputed values for variable \code{j}. If \code{y} is a \code{mids} object
#'then \code{imp[[j]]} equals \code{rbind(x$imp[[j]], y$imp[[j]])}; otherwise
#'the original data of \code{y} will be copied into this list, including the
#'missing values of \code{y} then \code{y} is not imputed.
#'Component \code{chainMean} is set to \code{NA} if \code{y} is \code{mids}, otherwise \code{x$chainMean} (which is
#' probably not reliable).
#'Component \code{chainVar} is set to \code{NA} if \code{y} is \code{mids}, otherwise \code{x$chainMean} (which is
#' probably not reliable).
#'Component \code{prepared} is set to \code{x$prepared} unless \code{y} is \code{mids}, in which case
#' the \code{prepared$data} elements are combined.
#'Component \code{loggedEvents} is set to \code{x$loggedEvents}.
#'Other elements are simply copied from \code{x}.
#'@note The results of imputations on subsets of the data, even when combined, will not be the same
#' as the result of imputation on the entire data.
#'@author Karin Groothuis-Oudshoorn, Stef van Buuren, 2009
#'@seealso \code{\link{cbind.mids}}, \code{\link{ibind}}, \code{\link[=mids-class]{mids}}
#'@references van Buuren S and Groothuis-Oudshoorn K (2011). \code{mice}:
#'Multivariate Imputation by Chained Equations in \code{R}. \emph{Journal of
#'Statistical Software}, \bold{45}(3), 1-67.
#'\url{http://www.jstatsoft.org/v45/i03/}
#'@keywords manip
#'@export
#
# TODO: ditch these semantics unless there is a compelling reason to keep them.
# I think it would be much more useful to require that y and ... all be mids objects,
# and to combine them assuming they represent independently imputed datasets.
# That is, they refer to different imputed datasets (associated with the m argument to mice)
# based on the original dataset, rather than representing runs on different subsets of the data.
# The current semantics when y is not a mids object don't seem very useful; I would eliminate them.
# RB
rbind.mids <- function(x, y, ...) {
    # This function combines x and y rowwise into a single midsobject.
    # x is a midsobject; y should be
    # 1) a dataframe with the same (number of) columns as x$data; in this case y
    # is combined with x$data with rbind() and the list elements of x are adjusted.
    # 2) y is a midsobject with the same underlying multiple imputation model as x but based on a
    # different data subset (with exactly the same variable(names) as in x$data). In this case the data is
    # combined with rbind and the other list elements of the midsobject are adjusted. Beware that
    # imputations in y are generated independently of x and by combining them could be dangerous.
    #
    # It is allowed to combine more than two objects when y is not a midsobject.
    # KO 08/09.
    call <- match.call()
    if (!is.mids(y))
        y <- rbind.data.frame(y, ...)

    # Then y is matrix, y is converted into a dataframe.
    if (is.matrix(y))
        y <- as.data.frame(y)

    varnames <- c(dimnames(x$data)[[2]])
    prepared <- x$prepared

    if (is.data.frame(y)) {
        if (ncol(y) != ncol(x$data))
            stop("The two datasets do not have the same number of columns\n")

        # The data in x (x$data) and y are combined together.
        data <- rbind(x$data, y)

        # count the number of missing data in y and add them to x$nmis.
        nmis <- x$nmis + colSums(is.na(y))

        # The original data of y will be copied into the multiple imputed dataset, including the missing values of y.
        r <- (!is.na(y))
        imp <- vector("list", ncol(y))
        for (j in visitSequence) {
            imp[[j]] <- rbind(x$imp[[j]], as.data.frame(matrix(NA, nrow = sum(!r[, j]), ncol = x$m, dimnames = list(row.names(y)[r[,
                                                                                                                                   j] == FALSE], 1:m))))
        }
        names(imp) <- varnames

        # Take chain stats from x, but that doesn't really seem right.
        # Also, the previous documentation claimed these values were always NA.
        # RB
        chainMean = x$chainMean
        chainVar = x$chainVar

        z <- list(data = data, nmis = nmis, imp = imp,
                  chainMean = chainMean, chainVar = chainVar, prepared = prepared)
    }

    if (is.mids(y)) {
        if (ncol(y$data) != ncol(x$data))
            stop("The two datasets do not have the same number of columns.\n")
        if (!all(c(dimnames(x$data)[[2]]) == c(dimnames(y$data)[[2]])))
            stop("The two datasets do not have the same variable(names).\n")
        if (!(x$m == y$m))
            stop("The number of imputations differ between x and y.\n")

        if (!identical(x$control == y$control))
            warning("Control values are not equal in x and y; will use values from x.\n")
        if (!all(x$method == y$method))
            warning("Methods vector is not equal in x and y; will use values from x.\n")
        if (!all(x$predictorMatrix == y$predictorMatrix))
            warning("Predictormatrix is not equal in x and y; will use values from x.\n")
        if (!all(x$visitSequence == y$visitSequence))
            warning("Visitsequence is not equal in x and y; will use values from x.\n")
        if (!all(x$form == y$form))
            warning("Formulae not equal in x and y; will use values from x.\n")
        if (!all(x$post == y$post))
            warning("The post vector is not equal in x and y; will use values from x.\n")


        # The data in x (x$data) and y are combined together.
        data <- rbind(x$data, y$data)
        prepared$data <- rbind(prepared$data, y$prepared$data)

        # count the number of missing data in y and add them to x$nmis.
        nmis <- x$nmis + y$nmis

        # The original data of y will be binded into the multiple imputed dataset, including the imputed values of y.
        imp <- vector("list", ncol(x$data))
        for (j in 1:ncol(x$data)) {
            imp[[j]] <- rbind(x$imp[[j]], y$imp[[j]])
        }
        names(imp) <- varnames

        # no meaningful chain statistics
        # chainMean an average of original chain means seems safe.  RB
        chainMean = NA
        chainVar = NA

        z <- list(data = data, m = m, nmis = nmis, imp = imp,
                  chainMean = chainMean, chainVar = chainVar, prepared = prepared)
    }

    # Call is a vector, with first argument the mice statement and second argument the call to cbind.mids.
    call <- match.call()
    z$call <- c(x$call, call)

    # Copy data from x
    z$m <- x$m
    z$iteration <- x$iteration
    z$control <- x$control
    z$method <- x$method
    z$form <- x$form
    z$extra <- x$extra
    z$post <- x$post
    z$predictorMatrix <- x$predictorMatrix
    z$visitSequence <- x$visitSequence
    z$seed <- x$seed
    z$lastSeedvalue <- x$lastSeedvalue
    z$iteration <- x$iteration

    # It would probably be better to combine the x and y loggedEvents if y is mids.  RB
    z$loggedEvents <- x$loggedEvents


    oldClass(z) <- "mids"
    return(z)
}
