#' Ordinal data regression using the Ordered Stereotype Model (OSM).
#'
#' Fit a regression model to an ordered factor response. The model is NOT a
#' logistic or probit model because the link function is not the logit, but the
#' link function is log-based.
#'
#' This function should be used in a very similar way to \code{MASS::polr}, and
#' some of the arguments are the same as \code{polr}, but the ordinal model used
#' here is less restrictive in its assumptions than the proportional odds model.
#' However, it is still parsimonious i.e. it uses only a small number of
#' additional parameters compared with the proportional odds model.
#'
#' @param formula a formula expression as for regression models, of the form
#'   response ~ predictors. The response should be a factor (preferably an
#'   ordered factor), which will be interpreted as an ordinal response, with
#'   levels ordered as in the factor. The model must have an intercept: attempts
#'   to remove one will lead to a warning and be ignored. An offset may be used.
#'   See the documentation of formula for other details.
#'
#' @param data a data frame, list or environment in which to interpret the
#'   variables occurring in \code{formula}.
#'
#' @param weights optional case weights in fitting. Default to 1.
#'
#' @param start initial values for the parameters. See the Details section for
#'   information about this argument.
#'
#' @param ... additional arguments to be passed to optim, most often a control
#'   argument.
#'
#' @param subset expression saying which subset of the rows of the data should
#'   be used in the fit. All observations are included by default.
#'
#' @param na.action a function to filter missing data.
#'
#' @param Hess logical for whether the Hessian (the observed information matrix)
#'   should be returned.
#'
#' @param model logical for whether the model matrix should be returned.
#'
#' @details This model is the \emph{ordered stereotype} model (Anderson 1984,
#' Agresti 2010)
#'
#' It is more flexible than the proportional odds model but only adds a handful
#' of additional parameters. It is not a cumulative model, being instead defined
#' in terms of the relationships between each of the higher categories and the
#' lowest category that is treated as the reference category.
#'
#' Each of the higher categories has its own intercept term, mu_k, which is
#' similar to the zeta parameters in \code{polr}, but in the OSM each higher
#' category also has its own scaling parameter, phi_k, which adjusts the effect
#' of the covariates on the response. This allows the effect of the covariates
#' on the response to be slightly different for each category of the response,
#' thus making the model more flexible than the proportional odds model.
#'
#' The final set of parameters are coefficients for each of the covariates, and
#' these are equivalent to the coefs in \code{polr}. Higher or more positive
#' values of the coefficients increases the probability of the response being in
#' the higher categories, and lower or more negative values of the coefficients
#' increase the probability of the response being in the lower categories.
#'
#' The overall model takes the following form:
#'
#' log(P(Y = k | X)/P(Y = 1 | X)) = mu_k + phi_k*beta_vec^T x_vec
#'
#' for k = 2, ..., q, where x_vec is the vector of covariates for the
#' observation Y.
#'
#' mu_1 is fixed at 0 for identifiability of the model, and the phi_k parameters
#' are constrained to be ordered (giving the model its name) in the following
#' way:
#'
#' 0 = phi_1 <= phi_2 <= ... <= phi_k <= ... <= phi_q = 1.
#'
#' (The unordered stereotype model restricts phi_1 and phi_q but allows the
#' remaining phi_k to be in any order, and this is suitable for fitting the
#' model for nominal data. However, this package does not provide that option,
#' as it is already available in other packages which can fit the stereotype
#' model.)
#'
#' After fitting the model, the estimated values of the intermediate phi_k
#' values indicate a suitable numerical spacing of the ordinal response
#' categories that is based on the data. The spacings indicate how much distinct
#' information each of the corresponding levels provide. For example, if you
#' have five response categories and the fitted phi values are \code{(0, 0.04,
#' 0.6, 0.62, 1)} then this indicates that levels 1 and 2 provide very similar
#' information about the effect of the covariates on the response, and levels 3
#' and 4 provide very similar information as each other. The meaning of this is
#' that you could simplify the response by combining levels 1 and 2 and
#' combining levels 3 and 4 (i.e. reduce the levels to 1, 3 and 5) and you would
#' still be able to estimate the beta coefficients with similar accuracy.
#'
#' Another use for the phi_k values is that if you want to carry out further
#' analysis of the response, treating it as a numerical variable, then the phi
#' values are a better choice of numerical values for the response categories
#' than the default values 1 to q.
#'
#' \strong{\code{start}} argument values: \code{start} is a vector of start
#' values for estimating the model parameters.
#'
#' The first part of the \code{start} vector is starting values for the
#' coefficients of the covariates, the second part is starting values for the mu
#' values (per-category intercepts), and the third part is starting values for
#' the raw parameters used to construct the phi values.
#'
#' The length of the vector is [number of covariate terms] + [number of
#' categories in response variable - 1] + [number of categories in response
#' variable - 2]. Every one of the values can take any real value.
#'
#' The second part is the starting values for the mu_k per-category intercept
#' parameters, and since mu_1 is fixed at 0 for identifiability, the number of
#' non-fixed mu_k parameters is one fewer than the number of categories.
#'
#' The third part of the starting vector is a re-parametrization used to
#' construct starting values for the estimated phi parameters such that the phi
#' parameters observe the ordering restriction of the ordered stereotype model,
#' but the raw parameters are not restricted which makes it easier to optimise
#' over them. phi_1 is always 0 and phi_q is always 1 (where q is the number of
#' response categories). If the raw parameters are u_2 up to u_(q-1), then phi_2
#' is constructed as expit(u_2), phi_3 is expit(u_2 + exp(u_3)), phi_4 is
#' expit(exp(u_3) + exp(u_4)) etc. which ensures that the phi_k values are
#' non-decreasing.
#'
#' This code was adapted from file MASS/R/polr.R
#' copyright (C) 1994-2013 W. N. Venables and B. D. Ripley
#' Use of transformed intercepts contributed by David Firth
#' The osm and osm.fit functions were written by Louise McMillan, 2020.
#'
#' This program is free software; you can redistribute it and/or modify it under
#' the terms of the GNU General Public License as published by the Free Software
#' Foundation; either version 2 or 3 of the License (at your option).
#'
#' This program is distributed in the hope that it will be useful, but WITHOUT
#' ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#' FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
#' details.
#'
#' A copy of the GNU General Public License is available at
#' http://www.r-project.org/Licenses/
#'
#'
#' @returns An object of class \code{"osm"}.  This has components
#'
#'   \code{beta} the coefficients of the covariates, with NO intercept.
#'
#'   \code{mu} the intercepts for the categories.
#'
#'   \code{phi} the score parameters for the categories (restricted to be
#'   ordered).
#'
#'   \code{deviance} the residual deviance.
#'
#'   \code{fitted.values} a matrix of fitted values, with a column for each
#'   level of the response.
#'
#'   \code{lev} the names of the response levels.
#'
#'   \code{terms} the \code{terms} structure describing the model.
#'
#'   \code{df.residual} the number of residual degrees of freedom, calculated
#'   using the weights.
#'
#'   \code{edf} the (effective) number of degrees of freedom used by the model
#'
#'   \code{n, nobs} the (effective) number of observations, calculated using the
#'   weights.
#'
#'   \code{call} the matched call.
#'
#'   \code{convergence} the convergence code returned by \code{optim}.
#'
#'   \code{niter} the number of function and gradient evaluations used by
#'   \code{optim}.
#'
#'   \code{eta}
#'
#'   \code{Hessian} (if \code{Hess} is true).  Note that this is a numerical
#'   approximation derived from the optimization proces.
#'
#'   \code{model} (if \code{model} is true), the model used in the fitting.
#'
#'   \code{na.action} the NA function used
#'
#'   \code{xlevels} factor levels from any categorical predictors
#'
#' @references Agresti, A. (2010). \emph{Analysis of ordinal categorical data} (Vol. 656). John Wiley & Sons.
#'
#' @references Anderson, J. A. (1984). Regression and ordered categorical variables. \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, 46(1), 1-22.
#'
#' @seealso [MASS::polr()]
#'
#' @importFrom stats .getXlevels binomial glm.fit model.matrix model.offset
#'   model.response model.weights deviance
#'
#' @export
osm <- function(formula, data, weights, start, ..., subset,
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
    num_beta <- ncol(x)

    ## Drop the (Intercept) column from the model matrix
    if(xint > 0L) {
        x <- x[, -xint, drop = FALSE]
        num_beta <- num_beta - 1L
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
            num_beta <- ncol(x)
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
        spacing <- logit((1L:(llev-1))/(llev)) # just a guess
        gammas <- -coefs[1L] + spacing - spacing[q1]

        ## Also generate starting values for phi, assuming equal spacing in
        ## the space of phi and converting using the logit link to the space
        ## of the auxiliary variable u
        startingphi <- (1:(llev-1))/(llev)
        u2 <- logit(startingphi)[1]
        us <- log(diff(logit(startingphi)))

        ## Construct the full starting values vector, using the fact that
        ## coefs[1L] has already been incorporated into the gammas object
        start <- c(coefs[-1L], gammas, u2, us)
    } else if(length(start) != num_beta + (llev-1) + (llev-2))
        stop("'start' is not of the correct length")

    ## Now run the fitting
    ans <- osm.fit(x, y, wt, start, offset, hessian = Hess, ...)

    ## Extract parts of the fitted model, to use when calculating the fitted
    ## values for each observation
    ## "res" is the output object from optim(), which contains the hessian
    ## object if requested
    beta <- ans$beta
    mu <- ans$mu
    phi <- ans$phi
    u <- ans$u
    res <- ans$res
    deviance <- ans$deviance
    logLik <- ans$logLik

    ## Calculate the fitted values of each observation, which are the probabilities
    ## of getting each of the levels of the response
    eta <- if(num_beta) offset + drop(x %*% beta) else offset + rep(0, n)
    fitted <- matrix(1,n,llev)
    for (k in 2:(llev)) {
        fitted[,k] <- exp(pmax(pmin(50,mu[k]+phi[k]*eta),-100))
    }
    fitted <- fitted/rowSums(fitted)
    dimnames(fitted) <- list(row.names(m), lev)

    ## Count the number of calls to the function, some of which were used
    ## to numerically calculate the gradient
    niter <- c(f.evals = res$counts[1L], g.evals = res$counts[2L])

    ## Construct the output object
    fit <- list(beta = beta, coefficients = beta, mu = mu, phi = phi, u = u,
                deviance = deviance, logLik = logLik,
                fitted.values = fitted, lev = lev, terms = Terms,
                df.residual = sum(wt) - num_beta - (llev-1) - (llev-2),
                edf = num_beta + (llev-1) + (llev-2), n = sum(wt),
                nobs = sum(wt), call = match.call(),
                convergence = res$convergence, niter = niter, eta = eta)
    if(Hess) {
        dn <- c(names(beta), names(ans$mu[-1L]), names(ans$u))
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
        mu <- c(0, coefficients[num_beta + ind_mu_k])

        ## u are the auxiliary parameters for phi, which are used because they
        ## are free to take any value between -Inf and Inf, and don't have to be
        ## ordered, whereas the phi values, from their construction, will be
        ## increasing, and the end values will be 0 and 1
        u <- coefficients[num_beta + num_mu_k + ind_phi_k]
        phi <- c(0 , expit(cumsum(c(u[1L], exp(u[-1L])))), 1)
        eta <- offset
        if (num_beta) eta <- eta + drop(x %*% coefficients[ind_beta])

        ## Construct the probabilities of getting each level of the response for
        ## each observation
        theta <- matrix(1, nrow=n, ncol=num_mu_k+1)
        for (k in 2:(num_mu_k+1)) {
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
    num_beta <- ncol(x)
    ind_beta <- seq_len(num_beta)

    ## Count the number of levels of y, and calculate q and q2 (there will be
    ## q independent values of the alpha parameters, and q2 = q-1 independent
    ## values of phi)
    lev <- levels(y)
    if(length(lev) <= 2L) stop("response must have 3 or more levels")
    y <- unclass(y)
    num_mu_k <- length(lev) - 1L
    ind_mu_k <- seq_len(num_mu_k)
    num_phi_k <- length(lev) - 2L
    ind_phi_k <- seq_len(num_phi_k)

    ## Run optim, and extract the results
    res <- optim(start, fmin, method="L-BFGS-B", ...)
    beta <- res$par[seq_len(num_beta)]
    mu <- c(0,res$par[num_beta + ind_mu_k])
    u <- res$par[num_beta + num_mu_k + ind_phi_k]
    phi <- c(0,expit(cumsum(c(u[1L], exp(u[-1L])))),1)
    deviance <- 2 * res$value
    logLik <- -res$value ## Negative because fmin is the NEGATIVE loglikelihood that gets MINIMIZED
    names(mu) <- paste0("mu",lev) ## Includes mu1 = 0 as first value
    names(phi) <- paste0("phi",lev) ## Includes phi1 = 0 as first value and phiq = 1 as last value
    names(u) <- paste0("phiAux",lev[-c(1L,length(lev))])
    if(num_beta) names(beta) <- colnames(x)
    list(beta = beta, mu = mu, phi = phi, u = u, deviance = deviance, logLik = logLik, res = res)
}

#' @importFrom stats naprint
#' @export
print.osm <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control=NULL)
    }
    if(length(x$beta)) {
        cat("\nCoefficients beta:\n")
        print(x$beta, ...)
    } else {
        cat("\nNo coefficients\n")
    }
    cat("\nIntercepts mu:\n")
    print(x$mu, ...)
    cat("\nScore parameters phi:\n")
    print(x$phi, ...)
    cat("\nResidual Deviance:", format(x$deviance, nsmall=2L), "\n")
    cat("AIC:", format(x$deviance + 2*x$edf, nsmall=2L), "\n")
    cat("BIC:", format(x$deviance + x$edf*log(x$n), nsmall=2L), "\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
    if(x$convergence > 0)
        cat("Warning: did not converge as iteration limit reached\n")
    invisible(x)
}

#' @export
logLik.osm <- function(object, ...)
    structure(object$logLik, df = object$edf,
              nobs = object[["nobs"]], class = "logLik")

#' @importFrom stats nobs
#' @export
nobs.osm <- function(object, ...) object[["nobs"]]

#' @importFrom stats vcov update
#' @importFrom utils tail
#' @export
vcov.osm <- function(object, ...)
{
    ## jacobian of phi differentiated with respect to u, for delta method for
    ## phi which are reparametrized from u so that u were able to be
    ## unconstrainedly optimized
    ## Note that this is the Jacobian for phi_2, ..., phi_(q-1) i.e. the
    ## INDEPENDENT values of phi ONLY, not phi_1=0 or phi_q=1
    jacobian <- function(u) {
        m <- length(u)
        eu <- exp(u)
        csum <- cumsum(c(u[1L], eu[-1L]))
        ecsum <- exp(-csum)
        denom <- (1+ecsum)^2
        mat <- matrix(0, m, m)
        for (i in 1L:m) {
            for (j in 1:i) {
                if (j == 1) mat[i,j] <- ecsum[i]/denom[i]
                else mat[i,j] <- ecsum[i]*eu[i]/denom[i]
            }
        }
        mat
    }

    if(is.null(object$Hessian)) {
        message("\nRe-fitting to get Hessian\n")
        utils::flush.console()
        object <- update(object, Hess = TRUE,
                         start = c(object$beta, object$mu[-1], object$u))
    }
    # vc <- MASS::ginv(object$Hessian)

    vc <- solve(object$Hessian)
    u <- object$u
    u.ind <- length(object$beta) + length(object$mu[-1]) + seq_along(u)
    J <- jacobian(u)
    A <- diag(tail(u.ind,1))

    ## Apply delta method to get correct SE for phi
    A[u.ind, u.ind] <- J
    V <- A %*% vc %*% t(A)

    cat("Variance-covariance matrix for beta coefficients, mu intercepts and the independent elements of phi.\n")

    ## Change the names of the phi elements of the matrix
    vcov_dimnames <- dimnames(object$Hessian)
    vcov_dimnames[[1]][u.ind] <- paste0("phi",2:(length(u)+1))
    vcov_dimnames[[2]][u.ind] <- paste0("phi",2:(length(u)+1))
    structure(V,  dimnames = vcov_dimnames)
}

#' @importFrom stats vcov pnorm
#' @export
summary.osm <- function(object, digits = max(3, .Options$digits - 3),
                        correlation = FALSE, ...)
{
    llev <- length(object$lev)
    cc <- c(coef(object), object$mu[-1L], object$phi[-c(1,llev)])
    pc <- length(coef(object))
    coef <- matrix(0, pc+llev-1+llev-2, 4L, dimnames=list(names(cc),
                                                          c("Value", "Std. Error", "t value", "p value")))
    coef[, 1L] <- cc
    vc <- vcov(object)
    coef[, 2L] <- sd <- sqrt(diag(vc))
    coef[, 3L] <- coef[, 1L]/coef[, 2L]
    coef[, 4L] <- pnorm(abs(coef[, 3L]), lower.tail = FALSE) * 2
    object$coefficients <- coef
    object$pc <- pc
    object$digits <- digits
    if(correlation)
        object$correlation <- (vc[1:pc,1:pc]/sd[1:pc])/rep(sd[1:pc], rep(pc, pc))
    class(object) <- "summary.osm"
    object
}

#' @importFrom stats naprint
#' @export
print.summary.osm <- function(x, digits = x$digits, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control=NULL)
    }
    coef <- format(round(x$coefficients, digits=digits))
    pc <- x$pc
    llev <- length(x$lev)
    if(pc > 0) {
        cat("\nCoefficients:\n")
        print(x$coefficients[seq_len(pc), , drop=FALSE], quote = FALSE,
              digits = digits, ...)
    } else {
        cat("\nNo coefficients\n")
    }
    cat("\nIntercepts mu:\n")
    print(coef[(pc+1L):(pc+llev-1), , drop=FALSE], quote = FALSE,
          digits = digits, ...)
    cat("\nScore parameters phi:\n")
    print(coef[(pc+llev):(pc+llev*2-3), , drop=FALSE], quote = FALSE,
          digits = digits, ...)
    cat("\nResidual Deviance:", format(x$deviance, nsmall=2L), "\n")
    cat("AIC:", format(x$deviance + 2*x$edf, nsmall=2L), "\n")
    cat("BIC:", format(x$deviance + x$edf*log(x$n), nsmall=2L), "\n")

    if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
    if(!is.null(correl <- x$correlation)) {
        cat("\nCorrelation of Coefficients:\n")
        ll <- lower.tri(correl)
        correl[ll] <- format(round(correl[ll], digits))
        correl[!ll] <- ""
        print(correl[-1L, -ncol(correl)], quote = FALSE, ...)
    }
    invisible(x)
}
