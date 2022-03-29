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
# ==============================================================================
# Roxygen below:
#
#' Logistic regression using ordered stereotype model (OSM)
#'
#' Fits a logistic regression model to an ordered factor response using the
#' Ordered Stereotype Models (OSM), using the reparametrization of Fernandez et
#' al (2016) and enforcing the non-decreasing constraint on \code{phi_k} values.
#'
#' See Anderson (1984) for the original definition of the ordered stereotype
#' model, and see Fernandez et al. (2016) for the reparametrization of the
#' \code{phi} values. Also see Agresti (2010) for more information about the
#' stereotype model and other ordinal models.
#'
#' The model assumes that the response is an ordered variable Y with q levels.
#'
#' This model is a baseline logit model, in which the distributions of levels
#' k = 2,...,q of the response Y are each defined separately relative to the
#' baseline. But the model is more parsimonious than the fully flexible baseline
#' logit model, because the effects of the linear predictor on the logits of
#' levels k = 2,...,q are proportional to each other.
#'
#' The model takes the form
#'
#' log(P(Y = k|x)/P(Y = 1|x)) = mu_k + phi_k*eta for k = 2,...,q,
#'
#' with \code{eta} being the linear predictor, i.e. \code{beta*x} where
#' \code{x} are the covariates and \code{beta} are the coefficients. There is no
#' \code{beta_0} term, because that is absorbed into the class boundary parameters
#' \code{mu_k}.
#'
#' This function fits the values of \code{mu_k}, \code{phi_k} and \code{beta}.
#'
#' \code{mu_1} = 0 for identifiability, so there are \code{q-1} independent
#' values of \code{mu_k}.
#'
#' The \code{phi} values are constrained to be non-decreasing, with
#'
#' phi_1 = 0 <= phi_2 <= ... <= phi_k = 1
#'
#' so there are \code{q-2} independent but non-decreasing values of
#' \code{phi_k}. This is the ordered form of the stereotype model.
#'
#' This is unlike the stereotype models fitted with the package
#' \url{https://cran.r-project.org/web/packages/ordinalgmifs/index.html}, which
#' do NOT enforce the ordering constraint on the \code{phi} values, and merely
#' enforce the 0 and 1 lower and upper limits for \code{phi}.
#'
#' If the user specifies the starting values for \code{beta}, \code{mu} and
#' \code{phi}, they can be specified in a vector, in that order. There should be
#' \code{length(x)} independent values of \code{beta},  \code{q-1} independent
#' values of \code{mu}, and \code{q-2} independent values of \code{phi}. BUT the
#' starting values for \code{phi} do not correspond directly to \code{phi},
#' because \code{phi} is restricted to being increasing and between 0 and 1, so
#' instead the starting values are treated as elements \code{u[2:q-1]} of a
#' vector \code{u} which can be between \code{-Inf} and \code{+Inf}, and then
#'
#'     \code{phi[2] <- expit(u[2])} and
#'
#'     \code{phi[k] <- expit(u[2] + sum(exp(u[3:k])))} for k between 3 and q-1
#'
#'     \code{(phi[1] = 0 and phi[q] = 1)}.
#'
#' The fitted values of \code{phi_k} and \code{u_k} are both included in the output.
#'
#' @param formula: a formula expression as for regression models, of the form
#'   \code{response ~ predictors}. The response should be a factor (preferably
#'   an ordered factor), which will be interpreted as an ordinal response, with
#'   levels ordered as in the factor. The model must have an intercept, which
#'   will then be removed because it is replaced by the fitted mu values;
#'   attempts by the user to remove the intercept in advance will lead to a
#'   warning. An offset may be used.
#'   See the documentation of \link[stats]{formula} for other details.
#'
#' @param data: an optional data frame, list or environment in which to
#'   interpret the variables occurring in \code{formula}.
#'
#' @param weights: optional case weights in fitting. Default to 1.
#'
#' @param start: initial values for the parameters. This is in the format
#'   \code{c(beta coefficients, mu, u)}: see the Values section.
#'
#' @param ...: additional arguments to be passed to \link[stats]{optim}, most
#'   often a \code{control} argument.
#'
#' @param subset: expression saying which subset of the rows of the data should
#'   be used in the fit. All observations are included by default.
#'
#' @param na.action: a function to filter missing data.
#'
#' @param Hess: logical for whether the Hessian (the observed information
#'   matrix) should be returned. Use this if you intend to call \code{summary}
#'   or \code{vcov} on the fit.
#'
#' @param model: logical for whether the model matrix should be returned.
#'
#' @param method: logistic or probit or (complementary) log-log or cauchit
#'   (corresponding to a Cauchy latent variable).
#'
#' @return
#' A object of class "osm". This has components:
#' \describe{
#'   \item{\code{beta}}{the coefficients of the covariates within the linear predictor, which has no intercept.}
#'   \item{\code{mu}}{the intercepts for the class boundaries.}
#'   \item{\code{phi}}{the class-specific scoring parameters that scale the effect of the linear predictor.}
#'   \item{\code{u}}{the raw scoring parameters fitted by \code{optim}. See Details.}
#'   \item{\code{deviance}}{the residual deviance.}
#'   \item{\code{fitted.values}}{a matrix, with a column for each level of the response.}
#'   \item{\code{lev}}{the names of the response levels.}
#'   \item{\code{terms}}{the \code{terms} structure describing the model.}
#'   \item{\code{df.residual}}{the number of residual degrees of freedoms, calculated using the weights.}
#'   \item{\code{edf}}{the (effective) number of degrees of freedom used by the model}
#'   \item{\code{n, nobs}}{the (effective) number of observations, calculated using the weights. (\code{nobs} is for use by \link[MASS]{stepAIC}).}
#'   \item{\code{call}}{the matched call.}
#'   \item{\code{convergence}}{the convergence code returned by \code{optim}.}
#'   \item{\code{niter}}{the number of function and gradient evaluations used by \code{optim}.}
#'   \item{\code{eta}}{the linear predictor (including any offset).}
#'   \item{\code{Hessian}}{(if \code{Hess} is true), the Hessian for the final fit.
#'   Note that this is a numerical approximation derived from the optimization
#'   process.}
#'   \item{\code{model}}{(if \code{model} is true), the model frame used.}
#'   \item{\code{na.action}}{the function used to filter missing data.}
#'   \item{\code{xlevels}}{list of factor levels from any factors or character vectors in the \code{terms} of the model.}
#' }
#'
#' @references
#' Fernandez, D., Arnold, R., & Pledger, S. (2016). Mixture-based clustering for the ordered stereotype model. *Computational Statistics & Data Analysis*, 93, 46-75.
#' Anderson, J. A. (1984). Regression and ordered categorical variables. *Journal of the Royal Statistical Society: Series B (Methodological)*, 46(1), 1-22.
#' Agresti, A. (2010). *Analysis of ordinal categorical data* (Vol. 656). John Wiley & Sons.
#'
#' @export
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