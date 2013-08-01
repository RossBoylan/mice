# ------------------------------sampler-------------------------------

sampler <- function(p, data, m, imp, r, visitSequence, fromto, printFlag, ...)
# The sampler controls the actual Gibbs sampling iteration scheme This function is called by mice and mice.mids
#
# Authors: S van Buuren, K Groothuis-Oudshoorn Copyright (c) 1999-2008 TNO Quality of Life
{
    ## set up array for convergence checking
    from <- fromto[1]
    to <- fromto[2]
    maxit <- to - from + 1
    # extra holds dynamic state of the variable specific imputation.
    # extra is a list of lists.  The outer list is indexed by imputed dataset,
    # and the inner list is indexed by imputed variable.  Entries in the inner list
    # have entries that are initially NULL.  Specific imputation methods will provide an appropriately typed object
    # that they pass back, which in turn are passed into the imputation method on the next Gibbs step.
    extra <- lapply(seq(m), function(i) vector("list", p$nvar))

    ## mminfo holds a list whose elements are either NULL
    ## or the return values from mmexpand; see that function below for details.
    ## mminfo[[j]] gives information on the model matrix for outcome j.
    ## mminfo is only used by 2 level predictors
    mminfo <- vector("list", p$nvar)

    if (maxit > 0)
        chainVar <- chainMean <- array(0, dim = c(length(visitSequence), maxit, m), dimnames = list(dimnames(data)[[2]][visitSequence],
            1:maxit, paste("Chain", 1:m))) else chainVar <- chainMean <- NULL

    ## THE ITERATION MAIN LOOP: GIBBS SAMPLER
    if (maxit < 1)
        iteration <- 0 else {
        if (printFlag)
            cat("\n iter imp variable")
        for (k in from:to) {
            #begin k loop : iteration loop (Gibbs sampler)
            iteration <- k
            for (i in 1:m) {
                #begin i loop    : repeated imputation loop (generated dataset)
                if (printFlag)
                  cat("\n ", iteration, " ", i)

                ## fill the data with the last set of imputations for the particular generated dataset
                # j ranges over the indices of variables to impute
                for (j in visitSequence) p$data[!r[, j], j] <- imp[[j]][, i]

                ## augment the data with the actual dummy variables
                for (j in setdiff(p$visitSequence, visitSequence)) {
                  cat.columns <- p$data[, p$categories[j, 4]]
                  p$data[, (j:(j + p$categories[p$categories[j, 4], 2] - 1))] <- matrix((model.matrix(~cat.columns - 1)[, -1]),
                    ncol = p$categories[p$categories[j, 4], 2], nrow = nrow(p$data))
                }

                ## iterate once over the variables of the augmented model

                for (j in p$visitSequence) {
                  theMethod <- p$method[j]
                  vname <- dimnames(p$data)[[2]][j]

                  ## store current state
                  oldstate <- get("state", pos = parent.frame())
                  newstate <- list(it = k, im = i, co = j, dep = vname, meth = theMethod, log = oldstate$log)
                  assign("state", newstate, pos = parent.frame(), inherits = TRUE)

                  if (printFlag & theMethod != "dummy")
                    cat(" ", vname)
                  if (theMethod != "" & (!is.passive(theMethod)) & theMethod != "dummy") {
                    # for a true imputation method
                    if (substring(tolower(theMethod), 1, 2) != "2l") {
                      # RJ: for an non-multilevel imputation method
                      # RB: formula-based  specification
                      x <- model.matrix(p$form[[j]], p$data)
                      y <- p$data[, j]
                      ry <- r[, j]
                      nam <- vname
                      if (k == 1)
                        check.df(x, y, ry, ...)  # added 31/10/2012, throw warning for n(obs) < p case
                      f <- paste("mice.impute", theMethod, sep = ".")
                      keep <- remove.lindep(x, y, ry, ...)
                      x <- x[, keep, drop = FALSE]
                      innerReturn <- do.call(f, args = list(y, ry, x, control=p$control[[j]], extra=extra[[i]][[j]], ...))
                    } else {
                      # for a multilevel imputation method
                      predictors <- p$predictorMatrix[j, ] != 0
                      # RB: formula-based specification
                      x <- model.matrix(p$form[[j]], p$data)
                      y <- p$data[, j]
                      ry <- r[, j]
                      minfo <- mminfo[[j]]
                      # lazy initialization: compute first time it needed
                      if (is.null(minfo))
                          minfo <- mminfo[[j]] <- mmexpand(p$predictorMatrix[j,], x, p$data, p$form[[j]])
                      type <- minfo$mmtype
                      nam <- vname
                      if (k == 1)
                        check.df(x, y, ry, ...)  # added 31/10/2012, throw warning for n(obs) < p case
                      f <- paste("mice.impute", tolower(theMethod), sep = ".")
                      keep <- remove.lindep(x, y, ry, ...)
                      x <- x[, keep, drop = FALSE]
                      type <- type[keep]
                      innerReturn <- do.call(f, args = list(y, ry, x, type, control=p$control[[j]], extra=extra[[i]][[j]], ...))
                    }
                    if (inherits(innerReturn, "innerReturn")) {
                        extra[[i]][[j]] <- innerReturn$extra
                        imp[[j]][,i] <- innerReturn$imp
                    } else {
                        imp[[j]][,i] <- innerReturn
                    }
                    p$data[!r[, j], j] <- imp[[j]][, i]
                  } else if (is.passive(theMethod)) {
                      ## TODO: It looks as if this can be folded into the formula processing
                    imp[[j]][, i] <- model.frame(as.formula(theMethod), p$data[!r[, j], ])  #RJ - FIXED passive imputation: as.formula()
                    p$data[!r[, j], j] <- imp[[j]][, i]
                  } else if (theMethod == "dummy") {
                      stop("Oh my! dummy method still in use")
                    ## FEH
                    cat.columns <- p$data[, p$categories[j, 4]]
                    p$data[, (j:(j + p$categories[p$categories[j, 4], 2] - 1))] <- matrix((model.matrix(~cat.columns - 1)[,
                      -1]), ncol = p$categories[p$categories[j, 4], 2], nrow = nrow(p$data))
                    remove("cat.columns")
                  }

                  ## optional post-processing
                  cmd <- p$post[j]  # SvB Aug 2009
                  if (cmd != "") {
                    eval(parse(text = cmd))
                    p$data[!r[, j], j] <- imp[[j]][, i]
                  }
                }  # end j loop
            }  # end i loop
            k2 <- k - from + 1
            for (j in 1:length(visitSequence)) {
                jj <- visitSequence[j]
                if (!is.factor(data[, jj])) {
                  chainVar[j, k2, ] <- apply(imp[[jj]], 2, var)
                  chainMean[j, k2, ] <- colMeans(as.matrix(imp[[jj]]))  ##pm 04/02
                }
                if (is.factor(data[, jj])) {
                  for (mm in 1:m) {
                    nc <- as.integer(factor(imp[[jj]][, mm], levels = levels(data[, jj])))
                    chainVar[j, k2, mm] <- var(nc)
                    chainMean[j, k2, mm] <- mean(nc)
                  }
                }
            }
        }  # end iteration loop
        if (printFlag)
            cat("\n")
    }
    return(list(iteration = maxit, imp = imp, chainMean = chainMean, chainVar = chainVar, extra=extra))
}


check.df <- function(x, y, ry, ...) {
    # if needed, writes the df warning message to the log
    df <- sum(ry) - ncol(x)
    mess <- paste("df set to 1. # observed cases:", sum(ry), " # predictors:", ncol(x))
    if (df < 1)
        updateLog(out = mess, frame = 3)
}

remove.lindep <- function(x, y, ry, eps = 1e-04, maxcor = 0.99, allow.na = FALSE, ...) {
    if (ncol(x) == 0)
        return(NULL)
    if (eps <= 0)
        stop("\n Argument 'eps' must be positive.")
    xobs <- x[ry, , drop = FALSE]
    if (allow.na) {
        if (sum(ry) == 0) {
            # escape for columns with only missing data  SvB 10/3/2011
            updateLog(out = "No observed cases, predictor removal skipped", frame = 3)
            return(rep(TRUE, ncol(x)))
        }
    }
    yobs <- as.numeric(y[ry])
    keep <- unlist(apply(xobs, 2, var) > eps)
    keep[is.na(keep)] <- FALSE
    keep <- keep & suppressWarnings((unlist(apply(xobs, 2, cor, yobs)) < maxcor))
    if (all(!keep))
        warning("All predictors are constant or have too high correlation.")
    k <- sum(keep)
    cx <- cor(xobs[, keep, drop = FALSE], use = "all.obs")
    eig <- eigen(cx, symmetric = TRUE)
    ncx <- cx
    while (eig$values[k]/eig$values[1] < eps) {
        j <- (1:k)[order(abs(eig$vectors[, k]), decreasing = TRUE)[1]]
        keep[keep][j] <- FALSE
        ncx <- cx[keep[keep], keep[keep], drop = FALSE]
        k <- k - 1
        eig <- eigen(ncx)
    }
    if (!all(keep)) {
        out <- paste(dimnames(x)[[2]][!keep], collapse = ", ")
        updateLog(out = out, frame = 3)
    }
    return(keep)
}


## make list of collinear variables to remove
find.collinear <- function(x, threshold = 0.999, ...) {
    nvar <- ncol(x)
    x <- data.matrix(x)
    r <- !is.na(x)
    nr <- apply(r, 2, sum, na.rm = TRUE)
    ord <- order(nr, decreasing = TRUE)
    xo <- x[, ord, drop = FALSE]  ## SvB 24mar2011
    varnames <- dimnames(xo)[[2]]
    z <- suppressWarnings(cor(xo, use = "pairwise.complete.obs"))
    hit <- outer(1:nvar, 1:nvar, "<") & (abs(z) >= threshold)
    out <- apply(hit, 2, any, na.rm = TRUE)
    return(varnames[out])
}


updateLog <- function(out = NULL, meth = NULL, frame = 2) {
    s <- get("state", parent.frame(frame))
    r <- get("loggedEvents", parent.frame(frame))

    rec <- data.frame(it = s$it, im = s$im, co = s$co, dep = s$dep, meth = ifelse(is.null(meth), s$meth, meth), out = ifelse(is.null(out),
        "", out))

    if (s$log)
        rec <- rbind(r, rec)
    s$log <- TRUE
    assign("state", s, pos = parent.frame(frame), inherits = TRUE)
    assign("loggedEvents", rec, pos = parent.frame(frame), inherits = TRUE)
    return()
}

# Expand type information on original data to match the model matrix
# The "imputand" is the variable to be imputed.
# predictorType = row from the predictorMatrix
# mm = model matrix of predictors of the imputand
# form = formula for the imputand
# data = original data
#
# returns a list
#   mmtype = predictorType expanded to match model matrix
#       mmtype[i] gives the predictor type of modelMatrix[,i]
#   mmcolnames = column names of the model matrix
#   iMMtoDF is such that data[,iMMtoDF[i]] contributed to mm[,i] (DF for data frame)
# Note that the return values are specific to this particular imputand.
mmexpand <- function(predictorType, mm, data, form) {
    mmi <- reverse.map(mm, form, data)
    return(list(mmtype=predictorType[mmi], mmcolnames = colnames(mm), iMMtoDF = mmi))
}

# return a vector v such that data[,v[i]] contributed to mm[,i]
# mm = model matrix created by
# form = formula
# data = data 
# If multiple columns of the original data contributed to a mm column,
# simply use the first such column.
reverse.map <- function(mm, form, data){
    tt <- terms(form, data=data)
    ttf <- attr(tt, "factors")
    mmi <- attr(mm, "assign")
    # this depends on assign using same order as columns of factors
    # entries in mmi that are 0 (the intercept) are silently dropped
    ttf2 <- ttf[,mmi]
    # take the first row that contributes
    r <- apply(ttf2, 2, function(is) rownames(ttf)[is > 0][1])
    v <- match(r, colnames(data))
    # assume if there is an intercept it will be at start
    if (0 %in% mmi)
        v <- c(NA, v)
    return(v)
}

