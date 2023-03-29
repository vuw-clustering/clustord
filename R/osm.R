# This code was adapted from
# file MASS/R/polr.R
# copyright (C) 1994-2013 W. N. Venables and B. D. Ripley
# Use of transformed intercepts contributed by David Firth
# The osm and osm.fit functions were written by Louise McMillan, 2020.
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
#' @importFrom stats .getXlevels binomial glm.fit model.matrix model.offset model.response model.weights deviance
osm <-
    function(formula, data, weights, start, ..., subset,
             na.action, Hess = FALSE, model = TRUE)
    {
        ## Create an object containing the original function call, which will
        ## later be used to obtain the model frame
        ## "expand.dots=FALSE" means those parts are left as "..." instead of
        ## converted to named arguments
        m <- match.call(expand.dots = FALSE)

        ## Convert the data from a matrix to a data frame if required
        if(is.matrix(eval.parent(m$data))) m$data <- as.data.frame(data)

        ## Delete any input arguments not needed to create the model frame,
        ## including the ... arguments
        m$start <- m$Hess <- m$model <- m$... <- NULL

        ## Keep the relevant inner parts of the function call, but change it to
        ## a call to "model.frame()"
        m[[1L]] <- quote(stats::model.frame)

        ## Evaluate this call to model.frame() in the parent frame i.e. where
        ## osm() was called from. This gets the model frame that contains only
        ## the data rows in the subset, and only the variables required by the
        ## formula
        m <- eval.parent(m)

        ## Also attach the Terms object of the formula to the model frame object
        Terms <- attr(m, "terms")

        ## Now create the model matrix, i.e. just the predictors, not the response,
        ## and with each categorical predictor changed to multiple dummy variables
        ## in the manner required by the "contrasts" option
        x <- model.matrix(Terms, m)

        ## Get the column index of the (Intercept) column that was created by
        ## the call to model.matrix()
        xint <- match("(Intercept)", colnames(x), nomatch = 0L)

        ## Count rows and columns in the model matrix
        n <- nrow(x)
        pc <- ncol(x)

        ## Drop the (Intercept) column from the model matrix
        if(xint > 0L) {
            x <- x[, -xint, drop = FALSE]
            pc <- pc - 1L
        } else warning("an intercept is needed and assumed")

        ## Fetch the weights, or generate ones
        wt <- model.weights(m)
        if(!length(wt)) wt <- rep(1, n)

        ## Fetch the offsets, or generate zeros
        offset <- model.offset(m)
        if(length(offset) <= 1L) offset <- rep(0, n)

        ## Fetch the response variable, and its levels, and count them
        y <- model.response(m)
        if(!is.factor(y)) stop("response must be a factor")
        lev <- levels(y); llev <- length(lev)
        if(llev <= 2L) stop("response must have 3 or more levels")
        y <- unclass(y)
        q <- llev - 1L

        ## Generate starting values for optimization
        if(missing(start)) {
            # try logistic/probit regression on 'middle' cut to find starting
            # values for the coefficients of the predictors
            # q1 is the level at, or just before, halfway through the levels of y
            q1 <- llev %/% 2L

            ## y1 is a binary variable with y1=0 if y <= q1 and y1=1 if y > q1
            y1 <- (y > q1)

            ## Construct a new model matrix and add an intercept column to it
            X <- cbind(Intercept = rep(1, n), x)

            ## Now attempt to fit logistic regression to the binary response y1
            fit <- glm.fit(X, y1, wt, family = binomial(), offset = offset)
            if(!fit$converged)
                stop("attempt to find suitable starting values failed")
            coefs <- fit$coefficients
            if(any(is.na(coefs))) {
                warning("design appears to be rank-deficient, so dropping some coefs")
                keep <- names(coefs)[!is.na(coefs)]
                coefs <- coefs[keep]
                x <- x[, keep[-1L], drop = FALSE]
                pc <- ncol(x)
            }

            ## The other parameters are labelled as alphas in Agresti's definition
            ## of the proportional odds model. They are the base probabilities
            ## for each level of the response variable, and must be strictly increasing.
            ## Generate them initially assuming they're evenly spaced across the
            ## range, convert to the linear predictor space using the logit link
            ## and adjust them to incorporate the fact that the logistic
            ## regression fitting produced an intercept term coefs[1L] for the
            ## q1 level
            logit <- function(p) log(p/(1 - p))
            spacing <- logit((1L:q)/(q+1L)) # just a guess
            gammas <- -coefs[1L] + spacing - spacing[q1]

            ## Also generate starting values for phi, assuming equal spacing in
            ## the space of phi and converting using the logit link to the space
            ## of the auxiliary variable u
            startingphi <- (1:(q-1))/q
            u2 <- logit(startingphi)[1]
            us <- log(diff(logit(startingphi)))

            ## Construct the full starting values vector, using the fact that
            ## coefs[1L] has already been incorporated into the gammas object
            start <- c(coefs[-1L], gammas, u2, us)
        } else if(length(start) != pc + q)
            stop("'start' is not of the correct length")

        ## Now run the fitting
        ans <- osm.fit(x, y, wt, start, offset, hessian = Hess, ...)

        ## Extract parts of the fitted model, to use when calculating the fitted
        ## values for each observation
        ## "res" is the output object from optim(), which contains the hessian
        ## object if requested
        beta <- ans$beta
        mu <- c(0,ans$mu)
        phi <- ans$phi
        res <- ans$res

        ## Calculate the fitted values of each observation, which are the probabilities
        ## of getting each of the levels of the response
        eta <- if(pc) offset + drop(x %*% beta) else offset + rep(0, n)
        fitted <- matrix(1,n,q+1)
        for (k in 2:(q+1)) {
            fitted[,k] <- exp(pmax(pmin(50,mu[k]+phi[k]*eta),-100))
        }
        fitted <- fitted/rowSums(fitted)
        dimnames(fitted) <- list(row.names(m), lev)

        ## Count the number of calls to the function, some of which were used
        ## to numerically calculate the gradient
        niter <- c(f.evals = res$counts[1L], g.evals = res$counts[2L])

        ## Construct the output object
        fit <- list(beta = beta, mu = mu, phi = phi, deviance = deviance,
                    fitted.values = fitted, lev = lev, terms = Terms,
                    df.residual = sum(wt) - pc - q, edf = pc + q, n = sum(wt),
                    nobs = sum(wt), call = match.call(),
                    convergence = res$convergence, niter = niter, eta = eta)
        if(Hess) {
            dn <- c(names(beta), names(ans$mu), names(ans$u))
            H <- res$hessian
            dimnames(H) <- list(dn, dn)
            fit$Hessian <- H
        }
        if(model) fit$model <- m
        fit$na.action <- attr(m, "na.action")
        fit$xlevels <- .getXlevels(Terms, m)
        class(fit) <- "osm"
        fit
    }

#' @importFrom stats optim
osm.fit <- function(x, y, wt, start, offset, ...)
{
    ## Set up the function call to use in optim(), which extracts the parameters
    ## from the parameter vector and calculates the negative of the log-likelihood
    fmin <- function(coefficients) {
        mu <- c(0, coefficients[pc + ind_q])

        ## u are the auxiliary parameters for phi, which are used because they
        ## are free to take any value between -Inf and Inf, and don't have to be
        ## ordered, whereas the phi values, from their construction, will be
        ## increasing, and the end values will be 0 and 1
        u <- coefficients[pc + q + ind_q2]
        phi <- c(0 , expit(cumsum(c(u[1L], exp(u[-1L])))), 1)
        eta <- offset
        if (pc) eta <- eta + drop(x %*% coefficients[ind_pc])

        ## Construct the probabilities of getting each level of the response for
        ## each observation
        theta <- matrix(1, nrow=n, ncol=q+1)
        for (k in 2:(q+1)) {
            theta[,k] <- exp(pmax(pmin(50,mu[k]+phi[k]*eta),-100))
        }
        theta <- theta/rowSums(theta)

        ## Now calculate the components of the likelihood for each observation
        pr <- vapply(1:n, function(i) theta[i,y[i]], 1)

        ## Construct the negative log-likelihood
        if (all(pr > 0)) -sum(wt * log(pr)) else Inf
    }

    ## Count the number of rows and columns in the model matrix, and the number
    ## of predictors (including dummy variables, not the original categorical variables)
    n <- nrow(x)
    pc <- ncol(x)
    ind_pc <- seq_len(pc)

    ## Count the number of levels of y, and calculate q and q2 (there will be
    ## q independent values of the alpha parameters, and q2 = q-1 independent
    ## values of phi)
    lev <- levels(y)
    if(length(lev) <= 2L) stop("response must have 3 or more levels")
    y <- unclass(y)
    q <- length(lev) - 1L
    ind_q <- seq_len(q)
    q2 <- length(lev) - 2L
    ind_q2 <- seq_len(q2)

    ## Run optim, and extract the results
    res <- optim(start, fmin, method="L-BFGS-B", ...)
    beta <- res$par[seq_len(pc)]
    mu <- res$par[pc + ind_q]
    u <- res$par[pc + q + ind_q2]
    phi <- c(0,expit(cumsum(c(u[1L], exp(u[-1L])))),1)
    deviance <- 2 * res$value
    names(mu) <- paste(lev[1L], lev[-1L], sep="|")
    names(phi) <- lev
    names(u) <- paste0("phiAux",lev[-c(1L,length(lev))])
    if(pc) names(beta) <- colnames(x)
    list(beta = beta, mu = mu, phi = phi, u=u, deviance = deviance, res = res)
}