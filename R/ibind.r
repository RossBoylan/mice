#'Combine imputations fitted to the same data
#'
#'This function combines two \code{mids} objects \code{x} and \code{y} into a
#'single \code{mids} object. The two \code{mids} objects should have the same
#'underlying multiple imputation model and should be fitted on exactly the same
#'dataset. If the number of imputations in \code{x} is \code{m(x)} and in
#'\code{y} is \code{m(y)} then the combination of both objects contains
#'\code{m(x)+m(y)} imputations.
#'
#'@param x A \code{mids} object.
#'@param y A \code{mids} object.
#'@return An S3 object of class \code{mids}
#'@note Although the resulting object has \code{$prepared}, it is just a copy of \code{x$prepared}, not
#' a meaningful combination.  Use the values in the \code{mids} object when possible.
#'@author Karin Groothuis-Oudshoorn, Stef van Buuren, 2009
#'@seealso \code{\link[=mids-class]{mids}}, \code{\link{rbind.mids}}, \code{\link{cbind.mids}}
#'@keywords manip
#'@export
ibind <- function(x, y) {
    ## This function combines two midsobjects x and y into a
    ## single midsobject. The two midsobjects should have the same underlying multiple imputation model and should
    ## have exactly the same dataset. If the number of imputations in x is m(x) and y is m(y) then the combination of both
    ## objects contains m(x)+m(y) imputations.
    ## KO 08/09.

    call <- match.call()
    call <- c(x$call, call)
    prepared <- x$prepared

    if (!is.mids(y) & !is.mids(x))
        stop("Both x and y should be a midsobject\n")
    if ((!all(is.na(x$data) == is.na(y$data))) | (!all(x$data[!is.na(x$data)] == y$data[!is.na(y$data)])))
        stop("The data in x and y and/or their reponsepattern do not completely agree")
    if (!all(x$nmis == y$nmis))
        stop("Number of missings is not equal in x and y.\n")
    if (!identical(x$control == y$control))
        stop("Control values are not equal in x and y.\n")
    if (!all(x$method == y$method))
        stop("Methods vector is not equal in x and y \n")
    if (!all(x$predictorMatrix == y$predictorMatrix))
        stop("Predictormatrix is not equal in x and y\n")
    if (!all(x$visitSequence == y$visitSequence))
        stop("Visitsequence is not equal in x and y \n")
    if (!all(x$form == y$form))
        stop("Formulae are not all equal in x and y.\n")
    if (!all(x$post == y$post))
        stop("The post vector is not equal in x and y \n")

    varnames <- c(dimnames(x$data)[[2]])
    visitSequence <- x$visitSequence
    imp <- vector("list", ncol(x$data))
    for (j in visitSequence) {
        imp[[j]] <- cbind(x$imp[[j]], y$imp[[j]])
    }

    m <- (x$m + y$m)
    iteration <- max(x$iteration, y$iteration)

    chainMean <- chainVar <- array(NA, dim = c(length(visitSequence), iteration, m),
                                   dimnames = list(names(visitSequence),
                                   1:iteration, paste("Chain", 1:m)))
    for (j in 1:x$m) {
        chainMean[, 1:x$iteration, j] <- x$chainMean[, , j]
        chainVar[, 1:x$iteration, j] <- x$chainVar[, , j]
    }
    for (j in 1:y$m) {
        chainMean[, 1:y$iteration, j + x$m] <- y$chainMean[, , j]
        chainVar[, 1:y$iteration, j + x$m] <- y$chainVar[, , j]
    }

    z <- list(call = call, data = x$data, m = m, nmis = x$nmis, imp = imp,
              method = x$method, predictorMatrix = x$predictorMatrix,
              visitSequence = visitSequence, form = x$form, control = x$control,
              post = x$post, seed = x$seed,
              iteration = iteration, lastSeedvalue = x$lastSeedvalue,
              extra = x$extra,
              chainMean = chainMean, chainVar = chainVar, prepared = prepared)

    oldClass(z) <- "mids"
    return(z)
}
